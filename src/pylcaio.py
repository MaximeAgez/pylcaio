"""
pyLCAIO hybridizes the ecoinvent v3.x and exiobase v3.x databases and exports the resulting hybrid database to a
brightway2 project.

author: Maxime Agez
e-mail: maxime.agez@polymtl.ca
"""

import numpy as np
import pandas as pd
import pkg_resources
import logging
import brightway2 as bw2
import pymrio
from tqdm import tqdm
import json
import wurst
import wurst.searching as ws
import uuid


class Hybridize:
    def __init__(self, bw2_project_name, ecoinvent_db_name, ecoinvent_version, path_to_exiobase):
        """
        Class hybridized ecoinvent with exiobase.

        :param bw2_project_name: [str] the name of your brightway2 project where you have an ecoinvent database stored
                                    that you wish to hybridize
        :param ecoinvent_db_name: [str] the name of the ecoinvent database you wih to hybridize
        :param ecoinvent_version: [str] the version of your ecoinvent database. Supported versions = ['3.9', '3.9.1',
                                    '3.10', '3.10.1']
        :param path_to_exiobase: [str] the path to the folder where you have stored exiobase files
        """

        # set up logging tool
        self.logger = logging.getLogger('pylcaio')
        self.logger.setLevel(logging.INFO)
        self.logger.handlers = []
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.propagate = False

        # set up brightway project
        if bw2_project_name not in bw2.projects:
            raise KeyError("The brightway project name passed does not match with any existing brightway project.")
        bw2.projects.set_current(bw2_project_name)

        # check ecoinvent database exists in project
        if ecoinvent_db_name not in bw2.databases:
            raise KeyError(
                "The ecoinvent database name passed does not match with any existing database within your brightway2 "
                "project.")
        self.db_name = ecoinvent_db_name
        if ecoinvent_version not in ['3.9', '3.9.1', '3.10', '3.10.1']:
            raise KeyError(
                "Your ecoinvent version is not supported by the current version of pyLCAIO. Supported versions are "
                "'3.9', '3.9.1', '3.10' and '3.10.1'")
        self.ei_version = ecoinvent_version

        self.logger.info("Getting technosphere matrix...")
        # need to create an LCA object to get access to the technosphere matrix later on
        self.lca = bw2.LCA({bw2.Database(self.db_name).random(): 1}, list(bw2.methods)[0])
        self.lca.lci()

        # load exiobase through pymrio
        self.logger.info("Loading exiobase...")
        self.io = pymrio.parse_exiobase3(path_to_exiobase)
        # exiobase is in millions euros -> convert to euros
        self.io.satellite.S[9:] /= 1000000
        self.io.satellite.unit[self.io.satellite.unit == 'M.EUR'] = 'euro'
        # get the reference year of used exiobase version
        self.reference_year_IO = int(self.io.meta.description[-4:])

        # load file with all filter and concordance information
        self.filter = pd.read_excel(pkg_resources.resource_filename(__name__, '/Data/mappings/filters.xlsx'), None)

        # concordance between ecoinvent countries/regions and exiobase countries/regions
        with open(pkg_resources.resource_filename(__name__, '/Data/mappings/geo_conc.json'), 'r') as f:
            self.concordance_geos = json.load(f)

        self.ei_wurst = None
        self.H = None
        self.G = None
        self.A_io_f = None
        self.A_io_f_uncorrected = None

        self.apply_wurst()

    def apply_wurst(self):
        """ wurst is used to fasten the identification of processes with brightway2"""

        self.logger.info("Load ecoinvent as wurst object to speed up data treatment...")
        self.ei_wurst = wurst.extract_brightway2_databases(self.db_name, add_identifiers=True)

    def get_uncorrected_upstream_cutoff_matrix(self):
        """
        Obtain the uncorrected upstream cutoff matrix.

        Brief methodology description:
        1. Obtain the concordance between ecoinvent processes and exiobase sectors and identify which processes of
        ecoinvent are hybridized (in the Excel files of the Data folder of pylcaio)
        2. Loop through the "to-hybridize" processes and apply concordances, both sectoral and geographical
        3. Through the concordance, copy the corresponding exiobase sectors and apply the price (corrected for inflation)
        to get the uncorrected upstream cutoffs
        """

        # build the concordance matrix (both sector and geo concordance)
        self.H = pd.DataFrame(0, index=pd.MultiIndex.from_product([self.io.get_regions().tolist(),
                                                                   self.io.get_sectors().tolist()]),
                              columns=self.filter['Hybridized processes'].code, dtype='float')

        # loop through the difference processes to-be-hybridized
        for col in tqdm(self.H.columns, leave=True):
            act = bw2.Database(self.db_name).get(col)
            # extract location and corresponding exiobase sector
            geo = act.as_dict()['location']
            sector = self.filter['Hybridized processes'].loc[
                self.filter['Hybridized processes'].code == act.as_dict()['code'], 'exiobase_sector'].iloc[0]
            # if ecoinvent process location is in exiobase regions (US -> US)
            if geo in list(self.H.index.levels[0]):
                self.H.loc[(geo, sector), col] = 1
            # if it needs some mapping
            elif geo in self.concordance_geos.keys():
                # it's a country, not a region (e.g., AR -> WL)
                if type(self.concordance_geos[geo]) == str:
                    self.H.loc[(self.concordance_geos[geo], sector), col] = 1
                # it's a region (e.g., RNA -> CA + US)
                else:
                    # then we need to do some weighted averages based on production values of the x vector of exiobase
                    self.H.loc[:, col] = (self.io.x.loc(axis=0)[self.concordance_geos[geo], sector] /
                                          self.io.x.loc(axis=0)[self.concordance_geos[geo], sector].sum()).reindex(
                        self.H.index).loc[:, 'indout'].fillna(0)
            # special case for the dynamic region of ecoinvent: RoW
            else:
                covered_geos_for_product = []
                # get all processes producing the reference product (wurst is way faster than bw2.search())
                filtering = ws.get_many(self.ei_wurst, ws.equals('reference product', act.as_dict()['reference product']))
                # extract the location of these processes
                for dataset in filtering:
                    if dataset['code'] in list(self.filter['Hybridized processes'].code):
                        covered_geos_for_product.append(dataset['location'])
                # only keep unique ones
                covered_geos_for_product = set(covered_geos_for_product)
                # remove RoW and GLO from set
                covered_geos_for_product = covered_geos_for_product - {'RoW'} - {'GLO'}
                # convert potential regions in countries
                covered_countries_for_product = [self.concordance_geos[i] for i in covered_geos_for_product]
                # convert potential list of lists as lists
                covered_countries_for_product = [x for xs in covered_countries_for_product for x in xs]
                # apply the weighted average for relevant countries to H
                self.H.loc[:, col] = (self.io.x.loc(axis=0)[[i for i in self.H.index.levels[0] if
                                                        i not in covered_countries_for_product], sector] /
                                      self.io.x.loc(axis=0)[[i for i in self.H.index.levels[0] if
                                                        i not in covered_countries_for_product], sector].sum()).reindex(
                    self.H.index).loc[:, 'indout'].fillna(0)

        # TODO get dynamic inflation numbers
        inflation = 1.25

        # upstream cutoff is basically the product of concordance (both sectoral and geographical) and price
        self.A_io_f_uncorrected = self.io.A.dot(self.H *
                                                self.filter['Hybridized processes'].set_index('code').price *
                                                inflation)

        # add a multiindex level to columns, being the name of the ecoinvent database
        self.H = pd.concat([self.H], keys=[self.db_name], axis=1)
        self.A_io_f_uncorrected = pd.concat([self.A_io_f_uncorrected], keys=[self.db_name], axis=1)
        # reindex to add non-hybridized processes as columns of zeros
        self.H = self.H.T.reindex(self.lca.activity_dict.keys()).fillna(0).T
        self.A_io_f_uncorrected = self.A_io_f_uncorrected.T.reindex(self.lca.activity_dict.keys()).fillna(0).T

    def correct_double_counting(self, method='STAM'):
        """
        Correcting the occurring double counting.

        Brief methodology description:
        1. First step is to add inputs of non-hybridized processes to the inventories of the hybridized processes. For
        example, adding the inputs of a coating step to the production of a steel plate overall. This transformation
        is given in the end by the dataframe only_hyb_with_non_hyb_inventory_converted.
        2. We decide which inputs are relevant or not by implementing thresholds below which inputs are not considered
        covered in the process' description
        3. We determine the filters correcting for double counting and applying them to the uncorrected upstream cutoff
        matrix.

        :param method: [str] The method for double counting correction can be either "binary" or "STAM".
                            Default option = "STAM"
        """

        # go to dense matrix
        A_ff = pd.DataFrame(self.lca.technosphere_matrix.todense(), self.lca.activity_dict.keys(),
                            self.lca.activity_dict.keys())
        # change notation (from LCA to IO notations)
        A_ff = pd.DataFrame(np.eye(len(A_ff)), A_ff.index, A_ff.columns) - A_ff

        # only keep technosphere matrix with non-hybridized processes
        A_not_hyb = A_ff.copy('deep')
        A_not_hyb.loc(axis=0)[:, self.filter['Hybridized processes'].code] = 0

        # only keep technosphere matrix with hybridized processes
        A_hyb = A_ff.copy('deep')
        not_hyb = (self.filter['Market processes'].code.tolist() +
                   self.filter['Internal and activities'].code.tolist() +
                   self.filter['Empty and aggregated processes'].code.tolist())
        A_hyb.loc(axis=0)[:, not_hyb] = 0

        # need to save RAM
        del A_ff

        # invert the technosphere matrix with only non-hybridized processes
        inverse = pd.DataFrame(np.linalg.inv(np.eye(len(A_not_hyb.values)) - A_not_hyb.values), A_not_hyb.index,
                               A_not_hyb.columns)

        # need to save RAM
        del A_not_hyb

        # multiply the inverse with A_hyb to include the inputs of non-hybridized processes in the description of
        # hybridized processes
        only_hyb_with_non_hyb_inventory_converted = A_hyb.dot(inverse)

        # need to save RAM
        del A_hyb

        # what constitutes a legit input in an LCA description? Does a 1e-11kg of plastic input /kg of product
        # warrant confidence that the LCA description is full and should not be completed? Nope. We set the threshold
        # of a relevant input to 1e-7, to account for cases with very limited quantities of metals used in coating.
        not_capital_goods = ws.get_many(self.ei_wurst, ws.exclude(ws.contains('unit', 'unit')))
        not_capital_goods_codes = []
        # extract the location of these processes
        for dataset in not_capital_goods:
            not_capital_goods_codes.append((dataset['database'], dataset['code']))
        # apply the threshold
        only_hyb_with_non_hyb_inventory_converted.loc[not_capital_goods_codes] = \
        only_hyb_with_non_hyb_inventory_converted.loc[not_capital_goods_codes].mask(
            only_hyb_with_non_hyb_inventory_converted.loc[not_capital_goods_codes] < 1e-7, 0)

        # for capital goods, quantities are always extremely small due to the allocation over the lifetime, threshold
        # is thus put to 1e-16 instead
        only_hyb_with_non_hyb_inventory_converted[only_hyb_with_non_hyb_inventory_converted < 1e-16] = 0

        if method == 'binary':
            # convert the ecoinvent inputs covered per process to exiobase sectors
            lambda_filter_matrix = self.H.dot(only_hyb_with_non_hyb_inventory_converted)
            # only keep inputs that are NOT covered by ecoinvent (i.e., equal to zero)
            lambda_filter_matrix = lambda_filter_matrix.mask(lambda_filter_matrix > 0)
            # changes the 0s to 1s
            lambda_filter_matrix[lambda_filter_matrix == 0] = 1
            # NaNs are the places where covered inputs were, fill with 0s
            lambda_filter_matrix = lambda_filter_matrix.fillna(0)

            # apply binary correction to uncorrected upstream cutoffs
            self.A_io_f = self.A_io_f_uncorrected.multiply(lambda_filter_matrix)

        elif method == 'STAM':
            # convert the ecoinvent inputs covered per process to exiobase sectors
            lambda_filter_matrix = self.H.dot(only_hyb_with_non_hyb_inventory_converted)
            # only keep inputs that are NOT covered by ecoinvent (i.e., equal to zero)
            lambda_filter_matrix = lambda_filter_matrix.mask(lambda_filter_matrix > 0)
            # changes the 0s to 1s
            lambda_filter_matrix[lambda_filter_matrix == 0] = 1
            # NaNs are the places where covered inputs were, fill with 0s
            lambda_filter_matrix = lambda_filter_matrix.fillna(0)

            # need to check if a missing input is legit or due to an omission
            self.STAM_table = pd.read_excel(pkg_resources.resource_filename(
                __name__, '/Data/double_counting_correction/STAM_table.xlsx'), index_col=0)
            with open(pkg_resources.resource_filename(
                    __name__, '/Data/double_counting_correction/STAM_categories.txt'), 'r') as file:
                self.io_categories = eval(file.read())

            # G contains the filter from the STAM table
            self.G = pd.DataFrame(0, index=self.io.get_sectors(), columns=self.STAM_table.columns)
            for col in self.G.columns:
                self.G.loc[self.io_categories[col], col] = 1
            # extend from 200 sectors to the 9800 sectors of exiobase
            self.G = pd.concat([self.G] * len(self.io.get_regions()), axis=0)
            # add first level of multiindex (regions)
            self.G.index = pd.MultiIndex.from_product([self.io.get_regions(), self.io.get_sectors()],
                                                      names=['region', 'sector'])
            # we consider that inputs of food in exiobase sectors are due to canteens and such -> exclude as we consider
            # employees would eat another way if canteen was not there -> canteen not required for the manufacture for
            # the commodity/service
            self.remove_canteen = pd.read_excel(pkg_resources.resource_filename(
                    __name__, '/Data/double_counting_correction/force_canteen_out.xlsx'), index_col=0)

            # gamma represents the inputs that are deemed missing on purpose based on STAM
            gamma_filter_matrix = self.G.dot((self.STAM_table.mul(self.remove_canteen)).dot(self.G.T.dot(self.H)))


            # phi represents the presence of inputs of similar functionality within the LCA description
            # presence of one of these inputs (e.g., diesel) forces the value zero to all other functionally similar
            # EEIO inputs (e.g., gasoline or kerosene)
            phi_filter_matrix = pd.DataFrame(1, index=self.G.index, columns=self.lca.activity_dict.keys())
            # convert the used inputs into used categories
            categories_used_by_processes = self.G.T.dot(self.H.dot(only_hyb_with_non_hyb_inventory_converted))
            # force to zero the following categories
            for category in ['Liquid Fuels', 'Solid Fuels', 'Gaseous Fuels', 'Electricity/heat', 'Transport']:
                phi_filter_matrix.loc[[i for i in phi_filter_matrix.index if i[1] in self.io_categories[category]],
                                      categories_used_by_processes.loc[category] != 0] = 0

            # apply various filters to remove the double counting from the uncorrected upstream cutoffs matrix
            self.A_io_f = phi_filter_matrix.multiply(
                gamma_filter_matrix.multiply(lambda_filter_matrix.multiply(self.A_io_f_uncorrected)))

        # save RAM
        del self.ei_wurst
        del self.A_io_f_uncorrected
        del self.H
        del self.G
        del self.lca

    def import_exiobase_in_bw2(self, culling_threshold=1e-7):
        """
        Function import exiobase3 within the brightway2 database
        :param culling_threshold: the threshold (in euros) below which the inputs will not be added to the IO description
                                  in brightway2. This ensures a smoother run of calculations as exiobase is  not sparse
                                  at all. And do we care if 1e-11 euros of whatever input for the production 1 euro of
                                  a given product is lost? Nope.
        :return:
        """

        # first, create the biosphere3 of exiobase3 environmental flows
        db_biosphere_exiobase_name = 'biosphere3_exiobase'

        exio3_biosphere = {}

        for env_stressor in self.io.satellite.unit.index:
            code = str(uuid.uuid4())

            if ' - air' in env_stressor:
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "emission",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('air',),
                    "name": env_stressor,
                    "code": code
                }
            elif ' - water' in env_stressor:
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "emission",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('water',),
                    "name": env_stressor,
                    "code": code
                }
            elif ' - soil' in env_stressor:
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "emission",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('soil',),
                    "name": env_stressor,
                    "code": code
                }
            elif ('Cropland' in env_stressor or 'Forest area' in env_stressor or
                  'Other land Use' in env_stressor or 'Permanent pastures' in env_stressor or
                  'Infrastructure land' in env_stressor):
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "natural resource",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('natural resource', 'land'),
                    "name": env_stressor,
                    "code": code
                }
            elif 'Extraction' in env_stressor and (
                    'Crop residues' in env_stressor or 'Fishery' in env_stressor or
                    'Fodder crops' in env_stressor or 'Forestry' in env_stressor or
                    'Grazing' in env_stressor or 'Primary Crops' in env_stressor):
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "natural resource",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('natural resource', 'biotic'),
                    "name": env_stressor,
                    "code": code
                }
            elif 'Extraction' in env_stressor and (
                    'Fossil Fuel' in env_stressor or 'Metal Ores' in env_stressor or
                    'Non-Metallic Minerals' in env_stressor):
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "natural resource",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('natural resource', 'in ground'),
                    "name": env_stressor,
                    "code": code
                }
            elif 'Water ' in env_stressor:
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "natural resource",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('natural resource', 'in water'),
                    "name": env_stressor,
                    "code": code
                }
            elif 'Energy ' in env_stressor:
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "inventory indicator",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('inventory indicator', 'resource use'),
                    "name": env_stressor,
                    "code": code
                }
            elif 'Emissions nec - waste - undef' == env_stressor:
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "inventory indicator",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('inventory indicator', 'waste'),
                    "name": env_stressor,
                    "code": code
                }
            elif ('axes' in env_stressor or 'wages' in env_stressor or 'Operating surplus' in env_stressor or
                  'Employment' in env_stressor):
                exio3_biosphere[(db_biosphere_exiobase_name, code)] = {
                    "type": "economic",
                    "unit": self.io.satellite.unit.loc[env_stressor, 'unit'],
                    "categories": ('economic', 'primary production factor'),
                    "name": env_stressor,
                    "code": code
                }

        # reformat in convenient dictionary for search of code from name of environmental flow
        exio3_biosphere_bw2_finding_codes = {v['name']: v['code'] for k, v in exio3_biosphere.items()}

        # then, create the brightway2 activities for exiobase3
        db_exiobase_name = 'exiobase'

        exio3_bw2 = {
            (db_exiobase_name, uuid.uuid4().hex): {
                "type": "process",
                "unit": "euro",
                "location": i[0],
                "name": i[1],
                "reference product": i[1],
                "database": db_exiobase_name,
                "exchanges": []}
            for i in self.io.A.columns
        }

        # add codes, can't add in dict comprehension because uuid.uuid4().hex is not fix
        for sector in exio3_bw2.keys():
            exio3_bw2[sector]['code'] = sector[1]

        # reformat in convenient dictionary for search of code from name and location
        exio3_bw2_finding_codes = {(v['name'], v['location']): v['code'] for k, v in exio3_bw2.items()}

        self.logger.info("Formatting technosphere inputs...")

        # add the production exchanges to the activities of exiobase3
        for sector in exio3_bw2:
            exio3_bw2[sector]['exchanges'].append({
                "amount": 1,
                "name": exio3_bw2[sector]['name'],
                "location": exio3_bw2[sector]['location'],
                "input": (db_exiobase_name, exio3_bw2[sector]['code']),
                "type": "production",
                "unit": "euro",
                "flow": str(uuid.uuid4())
            })

        # use stack to change from dataframe 9800x9800 to series 96040000x1 (faster that way)
        df = self.io.A.stack(future_stack=True).stack(future_stack=True)

        # exiobase is really not sparse with a bunch of very small values, so we cull those below the given threshold
        df = df[abs(df) > culling_threshold]

        # go to dict for speed
        df_dict = df.to_dict()

        # populate the technosphere exchanges within the activities of exiobase3
        for exc in tqdm(df_dict, leave=True):
            # identify purchasing sector
            purchaser_code = exio3_bw2_finding_codes[exc[2], exc[3]]
            # identify selling sector
            seller_code = exio3_bw2_finding_codes[exc[1], exc[0]]

            exio3_bw2[(db_exiobase_name, purchaser_code)]['exchanges'].append({
                "amount": df_dict[exc],
                "name": exc[1],
                "location": exc[0],
                "input": (db_exiobase_name, seller_code),
                "type": "technosphere",
                "unit": "euro",
                "flow": str(uuid.uuid4())
            })

        # use stack to change from dataframe 1113x9800 to series 10907400x1
        dff = self.io.satellite.S.stack(future_stack=True).stack(future_stack=True)
        # only keep non-zero values
        dff = dff[dff != 0]
        # go to dict for speed
        dff_dict = dff.to_dict()

        self.logger.info("Formatting biosphere inputs...")

        # populate the biosphere exchanges within the activities of exiobase3
        for exc in tqdm(dff_dict, leave=True):
            # identify emitting/extracting sector
            purchaser_code = exio3_bw2_finding_codes[exc[1], exc[2]]
            # identify environmental stressor emitted/extracted
            env_stressor_code = exio3_biosphere_bw2_finding_codes[exc[0]]

            exio3_bw2[(db_exiobase_name, purchaser_code)]['exchanges'].append({
                "amount": dff_dict[exc],
                "name": exc[0],
                "input": (db_biosphere_exiobase_name, env_stressor_code),
                "type": 'biosphere',
                "unit": self.io.satellite.unit.loc[exc[0], 'unit'],
                "flow": str(uuid.uuid4())
            })

        # write the biosphere database
        self.logger.info("Writing biosphere3 database for exiobase...")
        bw2.Database(db_biosphere_exiobase_name).write(exio3_biosphere)

        # write the exiobase3 database
        self.logger.info("Writing exiobase database...")
        bw2.Database(db_exiobase_name).write(exio3_bw2)

        # save some RAM
        del self.io

    def import_hybridized_database(self):
        """Function imports the hybridized version of ecoinvent into brightway2."""

        self.logger.info("Copying existing ecoinvent database...")
        hybrid_ei_data = bw2.Database(self.db_name).load()

        # get codes of exiobase in dictionary for speed
        db_exiobase_name = 'exiobase'
        exio_find_code_easily = {(i.as_dict()['location'], i.as_dict()['name']): i.as_dict()['code'] for i in
                                 bw2.Database(db_exiobase_name)}

        self.logger.info("Reformat the upstream cutoffs matrix...")
        # change to pd.Series and remove null inputs to speed things up
        upstream_cutoffs = self.A_io_f.stack(future_stack=True).stack(future_stack=True)
        upstream_cutoffs = upstream_cutoffs[upstream_cutoffs != 0]

        self.logger.info("Remove irrelevant upstream cutoffs inputs...")
        # store prices in a df to use them to select irrelevant inputs
        prices = self.filter['Hybridized processes'].set_index('code').loc[:, 'price']

        # we want to remove some of the inputs because it takes a long time to write the database otherwise
        # so if an input is smaller than 0.001% of the total price of the product -> irrelevant input
        irrelevant_inputs = []
        for exc in tqdm(list(upstream_cutoffs.keys()), leave=True):
            if upstream_cutoffs[exc] < (prices.loc[exc[2]] * 0.00001):
                irrelevant_inputs.append(exc)
        # drop the irrelevant inputs
        upstream_cutoffs = upstream_cutoffs.drop(irrelevant_inputs)

        # change to dict format to speed things up
        upstream_cutoffs = upstream_cutoffs.to_dict()

        # save RAM
        del self.A_io_f

        self.logger.info("Add relevant upstream cutoffs inputs to ecoinvent processes...")
        # put the upstream cutoffs in the ecoinvent descriptions
        for exc in tqdm(list(upstream_cutoffs.keys()), leave=True):
            # get the key of the ecoinvent process to hybridize
            receiving_process = (exc[3], exc[2])
            # get the key of the upstream cutoff to add to the ecoinvent process
            sectoral_input_bought = exio_find_code_easily[(exc[0], exc[1])]
            # add the exchange to the dictionary
            hybrid_ei_data[receiving_process]['exchanges'].append({
                "amount": upstream_cutoffs[exc],
                "name": exc[1],
                "location": exc[0],
                "type": "technosphere",
                "unit": "euro",
                "flow": str(uuid.uuid4()),
                "input": (db_exiobase_name, sectoral_input_bought)
            })

            # Remove the processed key to free memory
            del upstream_cutoffs[exc]

        self.logger.info("Change name of the ecoinvent database to hybrid ecoinvent database...")
        # change the keys of processes from "ecoinvent" to "hybrid ecoinvent"
        hybrid_ei_data = {('hybrid ' + self.db_name, k[1]): v for k, v in hybrid_ei_data.items()}
        # change the keys of inputs and outputs from "ecoinvent" to "hybrid ecoinvent"
        for k in tqdm(hybrid_ei_data.keys(), leave=True):
            hybrid_ei_data[k]['database'] = 'hybrid ' + self.db_name
            for exc in hybrid_ei_data[k]['exchanges']:
                if exc['type'] in ['technosphere', 'production']:
                    if exc['input'][0] == self.db_name:
                        exc['input'] = ('hybrid ' + self.db_name, exc['input'][1])
                    # for exiobase database, no output key was explicitly defined
                    if 'output' in exc.keys():
                        if exc['output'][0] == self.db_name:
                            exc['output'] = ('hybrid ' + self.db_name, exc['output'][1])

        self.logger.info("Writing the hybrid database to brightway...")
        db = bw2.Database('hybrid ' + self.db_name)
        db.write(hybrid_ei_data)

    def import_lcia_method_for_hybrid_lca(self):
        """Since hybrid LCA mixes elementary flows from ecoinvent and exiobase, we need an LCIA method that covres
        both at the same time for convenience. That's what this function creates."""

        # loads IW+ for ecoinvent
        if '3.9' in self.ei_version:
            iw_ei = bw2.BW2Package.load_file(pkg_resources.resource_filename(__name__,
                '/Data/LCIA/impact_world_plus_21_brightway2_expert_version_ei39.6cd1745d7173fc689a3cc8c44fd3e41d.bw2package'))
        elif '3.10' in self.ei_version:
            iw_ei = bw2.BW2Package.load_file(pkg_resources.resource_filename(__name__,
                '/Data/LCIA/impact_world_plus_21_brightway2_expert_version_ei310.5535d12bedce3770ffef004e84229fd1.bw2package'))

        # loads IW+ for exiobase
        C_exio = pd.read_excel(pkg_resources.resource_filename(__name__,
                               '/Data/LCIA/impact_world_plus_2.1_expert_version_exiobase.xlsx'), index_col=0)
        # change to series (more efficient)
        C_exio = C_exio.stack()
        # remove null CFs
        C_exio = C_exio[C_exio != 0]
        # get the Total human health impact
        total_hh = C_exio.loc[[i for i in C_exio.index if 'DALY' in i[0]]].groupby(axis=0, level=1).sum()
        total_hh = pd.concat([total_hh], keys=['Total human health (DALY)'])
        # get the Total ecosystem quality impact
        total_eq = C_exio.loc[[i for i in C_exio.index if 'PDF.m2.yr' in i[0]]].groupby(axis=0, level=1).sum()
        total_eq = pd.concat([total_eq], keys=['Total ecosystem quality (PDF.m2.yr)'])
        # concat the totals with the original CFs
        C_exio = pd.concat([C_exio, total_hh, total_eq])
        # a dictionary to rapidly find the keys for exiobase elementary flows
        exio3_biosphere_bw2_finding_codes = {act.as_dict()['name']: act.as_dict()['code'] for act in
                                             bw2.Database('biosphere3_exiobase')}
        # write the method combining CFs for ecoivnent and exiobase
        for method in iw_ei:
            method['name'] = ('for hybrid ecoinvent'.join(method['name'][0].split('for ecoinvent')),
                              method['name'][1], method['name'][2])
            try:
                for emission in C_exio.loc[method['name'][2] + ' (' + method['metadata']['unit'] + ')'].index:
                    method['data'].append((('biosphere3_exiobase', exio3_biosphere_bw2_finding_codes[emission]),
                                           C_exio.loc[method['name'][2] + ' (' + method['metadata']['unit'] + ')'].loc[
                                               emission]))
            # KeyError is the category has no CFs in exiobase (e.g., Ionizing radiations -> no ionizing radiation
            # elementary flows in exiobase)
            except KeyError:
                pass
            # write in brightway2
            method_to_import = bw2.Method(method['name'])
            method_to_import.register()
            method_to_import.metadata['unit'] = method['metadata']['unit']
            method_to_import.write(method['data'])