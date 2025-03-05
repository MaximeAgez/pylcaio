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
import gzip
import pickle


class Hybridize_ecoinvent:
    def __init__(self, bw2_project_name, ecoinvent_db_name, ecoinvent_version, path_to_exiobase):
        """
        Class hybridizes ecoinvent with exiobase.

        :param bw2_project_name: [str] the name of your brightway2 project where you have an ecoinvent database stored
                                    that you wish to hybridize
        :param ecoinvent_db_name: [str] the name of the ecoinvent database you wish to hybridize
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
        if ecoinvent_version in ['3.9', '3.9.1']:
            self.ei_version = '3.9'
        elif ecoinvent_version in ['3.10', '3.10.1']:
            self.ei_version = '3.10'

        self.logger.info("Getting technosphere matrix...")
        # need to create an LCA object to get access to the technosphere matrix later on
        self.lca = bw2.LCA({bw2.Database(self.db_name).random(): 1}, list(bw2.methods)[0])
        self.lca.lci()

        # load exiobase through pymrio
        self.logger.info("Loading exiobase...")
        io = pymrio.parse_exiobase3(path_to_exiobase)
        # we do not store the whole io object as many of the attributes are useless to us and we are short on RAM
        self.x_io = io.x.copy()
        self.A_io = io.A.copy()
        self.sectors_io = io.get_sectors().tolist()
        self.regions_io = io.get_regions().tolist()
        self.io_units = io.satellite.unit.copy()
        self.io_units[self.io_units == 'M.EUR'] = 'euro'
        # get the reference year of used exiobase version
        self.reference_year_IO = int(io.meta.description[-4:])
        # save RAM
        del io

        # load capital matrix of exiobase (capital endogenization)
        with gzip.open(pkg_resources.resource_filename(
                __name__, '/Data/capitals_exiobase/K_cfc_pxp_exio3.8.2_'+str(self.reference_year_IO)+'.gz.pickle'),
                'rb') as file:
            K = pickle.load(file)
        # add capital matrix to technology matrix. No need to adjust satellite accounts because they're not used
        self.A_io += K

        # load file with all filter and concordance information
        self.filter = pd.read_excel(pkg_resources.resource_filename(
            __name__, '/Data/ecoinvent/ei'+self.ei_version+'/mappings/filters.xlsx'), None)

        # concordance between ecoinvent countries/regions and exiobase countries/regions
        with open(pkg_resources.resource_filename(
                __name__, '/Data/ecoinvent/ei'+self.ei_version+'/mappings/geo_conc.json'), 'r') as f:
            self.concordance_geos = json.load(f)

        self.covered_geos = None
        self.not_capital_goods_codes = None
        self.codes_to_names = None
        self.A_io_f = None
        self.A_io_f_uncorrected = None

        self.get_relevant_info()

    def get_relevant_info(self):
        """ wurst is used to fasten the identification of processes with brightway2"""

        self.logger.info("Formatting relevant data of ecoinvent for quick access...")

        # identify which geographies are covered for which product (useful for determining which countries are in RoW)
        self.covered_geos = {act.as_dict()['reference product']: [] for act in bw2.Database(self.db_name)}

        for act in bw2.Database(self.db_name):
            if (act.as_dict()['location'] not in self.covered_geos[act.as_dict()['reference product']]
                    and act.as_dict()['code'] in self.filter['Hybridized processes']):
                self.covered_geos[act.as_dict()['reference product']].append(act.as_dict()['location'])

        # identify which products in ecoinvent are not capital goods, useful for removing irrelevant inputs later
        self.not_capital_goods_codes = [act.as_dict()['code'] for act in bw2.Database(self.db_name) if
                                        act.as_dict()['unit'] != 'unit']

        # get dictionary to quickly get a product name from the process code
        self.codes_to_names = {act.as_dict()['code']: act.as_dict()['reference product'] for act in
                               bw2.Database(self.db_name)}

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

        self.logger.info("Getting the uncorrected upstream cutoff matrix...")

        # build the concordance matrix (both sector and geo concordance)
        concordance_matrix = pd.DataFrame(0, index=pd.MultiIndex.from_product([self.regions_io, self.sectors_io]),
                                          columns=self.filter['Hybridized processes'].code, dtype='float')

        # loop through the difference processes to-be-hybridized
        for col in tqdm(concordance_matrix.columns, leave=True):
            act = bw2.Database(self.db_name).get(col)
            # extract location and corresponding exiobase sector
            geo = act.as_dict()['location']
            sector = self.filter['Hybridized processes'].loc[
                self.filter['Hybridized processes'].code == act.as_dict()['code'], 'exiobase_sector'].iloc[0]
            # if ecoinvent process location is in exiobase regions (US -> US)
            if geo in list(concordance_matrix.index.levels[0]):
                concordance_matrix.loc[(geo, sector), col] = 1
            # if it needs some mapping
            elif geo in self.concordance_geos.keys():
                # it's a country, not a region (e.g., AR -> WL)
                if type(self.concordance_geos[geo]) == str:
                    concordance_matrix.loc[(self.concordance_geos[geo], sector), col] = 1
                # it's a region (e.g., RNA -> CA + US)
                else:
                    # then we need to do some weighted averages based on production values of the x vector of exiobase
                    concordance_matrix.loc[:, col] = (self.x_io.loc(axis=0)[self.concordance_geos[geo], sector] /
                                                      self.x_io.loc(axis=0)[self.concordance_geos[geo], sector].sum()
                                                      ).reindex(concordance_matrix.index).loc[:, 'indout'].fillna(0)
            # special case for the dynamic region of ecoinvent: RoW
            else:
                # get covered geographies (both countries and regions)
                covered_geos_for_product = set(self.covered_geos[self.codes_to_names[col]]) - {'RoW'} - {'GLO'}
                # convert potential regions in countries (RER = FR + DE + BE + ...)
                covered_countries_for_product = [self.concordance_geos[i] for i in covered_geos_for_product]
                # convert potential list of lists as lists
                covered_countries_for_product = set([x if type(xs) == list else xs for xs in
                                                     covered_countries_for_product for x in xs])
                # apply the weighted average for relevant countries to H
                concordance_matrix.loc[:, col] = (self.x_io.loc(axis=0)[
                                                      [i for i in concordance_matrix.index.levels[0] if
                                                       i not in covered_countries_for_product], sector] /
                                                  self.x_io.loc(axis=0)[
                                                      [i for i in concordance_matrix.index.levels[0] if
                                                       i not in covered_countries_for_product], sector].sum()).reindex(
                    concordance_matrix.index).loc[:, 'indout'].fillna(0)

        # TODO get dynamic inflation numbers
        inflation = 1.25

        # upstream cutoff is basically the product of concordance (both sectoral and geographical) and price
        self.A_io_f_uncorrected = self.A_io.dot(concordance_matrix *
                                                self.filter['Hybridized processes'].set_index('code').price *
                                                inflation)

        # add a multiindex level to columns, being the name of the ecoinvent database
        self.A_io_f_uncorrected.columns = pd.MultiIndex.from_arrays(
            [[self.db_name] * len(self.A_io_f_uncorrected.columns), self.A_io_f_uncorrected.columns])

        # save RAM
        del self.x_io

    def determine_covered_inputs(self):
        """
        To correct for double counting te first step is to determine which inputs are already covered by the LCA
        database. To do so, we need to propagate the inputs of non-hybridized processes into inventories of hybridized
        inventories. So that if the production of timber (hybridized) requires power sawing (non-hybridized) the power
        saw is ultimately associated in the production of timber inventory so that an input of Machinery is not added
        through hybridization.
        Also, once all inputs have been propagated to hybridized processes, we only keep relevant inputs. What do we
        mean by that? Well if you have a random 1e-11 kg input of steel somewhere in your LCA description is very unlikely
        to mean that steel was properly covered by LCA. So how do we determine those relevant inputs? Well we fixed them
        to 1e-7 to acccount for cases where really small quantities are actually relevant, i.e., the use of metals for
        coating purposes for example.
        For capital goods specifically, the threshold has been set lower to acccount for the small quantities in LCA
        descriptions due to the distribution over the long lifetime of capital goods. The thrshold is therefore 1e-16
        for capital goods.
        """

        self.logger.info("Determining the covered inputs for each process of the LCI database...")

        # we need an LCA object to get to the technosphere matrix of brightway2
        lca = bw2.LCA({bw2.Database(self.db_name).random(): 1}, list(bw2.methods)[0])
        lca.lci()

        # go to dense matrix
        A_ff = pd.DataFrame(lca.technosphere_matrix.todense(), lca.activity_dict.keys(),
                            lca.activity_dict.keys())
        # save RAM
        del lca
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
        self.covered_inputs = A_hyb.dot(inverse)

        # need to save RAM
        del A_hyb
        del inverse

        # remove irrelevant inputs for non capital goods
        self.covered_inputs.loc(axis=0)[:, self.not_capital_goods_codes] = (
            self.covered_inputs.loc(axis=0)[:, self.not_capital_goods_codes].mask(
                self.covered_inputs.loc(axis=0)[:, self.not_capital_goods_codes] < 1e-7, 0))

        # remove irrelevant inputs for capital goods
        self.covered_inputs[self.covered_inputs < 1e-16] = 0

        # remove first level multiindex
        self.covered_inputs = self.covered_inputs.droplevel(0)

        # only keep non-null ecoinvent inputs (propagated ones became null)
        self.covered_inputs = self.covered_inputs.loc[self.covered_inputs.sum(1)[self.covered_inputs.sum(1) != 0].index]

        # convert ecoinvent inputs into exiobase sectoral inputs
        self.covered_inputs.index = [
            self.filter['Hybridized processes'].set_index('code').loc[i, 'exiobase_sector'] for i in
            self.covered_inputs.index]

        # groupby the sectoral inputs
        self.covered_inputs = self.covered_inputs.groupby(self.covered_inputs.index).sum()

        # only keep hybridized processes as columns
        self.covered_inputs = self.covered_inputs.loc(axis=1)[:, self.filter['Hybridized processes'].code]

        # convert codes to names
        self.covered_inputs.columns = [bw2.Database(i[0]).get(i[1]).as_dict()['reference product'] for i in
                                       self.covered_inputs.columns]

        # groupby columns
        self.covered_inputs = self.covered_inputs.groupby(axis=1, level=0).sum()

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

        # change size of covered_inputs to match exiobase sectors
        multi_index = pd.MultiIndex.from_product([self.regions_io, self.covered_inputs.index])
        self.covered_inputs = pd.concat([self.covered_inputs] * len(self.regions_io), axis=0)
        self.covered_inputs.index = multi_index
        self.covered_inputs = self.covered_inputs.reindex(self.A_io.index).fillna(0)

        # only keep inputs that are NOT covered by ecoinvent (i.e., equal to zero)
        self.covered_inputs = self.covered_inputs.mask(self.covered_inputs > 0)
        # changes the 0s to 1s
        self.covered_inputs[self.covered_inputs == 0] = 1
        # NaNs are the places where covered inputs were, fill with 0s
        self.covered_inputs = self.covered_inputs.fillna(0)

        # get the lambda filter matrix for regioinvent
        lambda_filter_matrix = pd.DataFrame(0, self.covered_inputs.index, self.A_io_f_uncorrected.columns, dtype=float)
        hybridized = self.filter['Hybridized processes'].code.tolist()
        for col in tqdm(lambda_filter_matrix.columns, leave=True):
            if col[1] in hybridized:
                lambda_filter_matrix.loc[:, col] = self.covered_inputs.loc[:, self.codes_to_names[col[1]]]

        # save RAM
        del self.covered_inputs
        del hybridized

        if method == 'binary':
            self.logger.info("Applying the various double counting correction matrices to the uncorrected upstream "
                             "cutoff matrix...")
            self.A_io_f = self.A_io_f_uncorrected.multiply(lambda_filter_matrix)

        elif method == 'STAM':

            self.logger.info("Getting the gamma matrix for double counting correction...")

            # need to check if a missing input is legit or due to an omission
            STAM_table = pd.read_excel(pkg_resources.resource_filename(
                __name__, '/Data/ecoinvent/ei'+self.ei_version+'/double_counting_correction/STAM_table.xlsx'),
                index_col=0)
            with open(pkg_resources.resource_filename(
                    __name__, '/Data/ecoinvent/ei'+self.ei_version+'/double_counting_correction/STAM_categories.txt'),
                    'r') as file:
                io_categories = eval(file.read())

            # STAM_df converts STAM_table to a format with exiobase sectors as indexes
            STAM_df = pd.DataFrame(0, index=self.sectors_io, columns=STAM_table.columns)
            for col in STAM_df.columns:
                STAM_df.loc[io_categories[col], col] = 1
            # extend from 200 sectors to the 9800 sectors of exiobase
            STAM_df = pd.concat([STAM_df] * len(self.regions_io), axis=0)
            # add first level of multiindex (regions)
            STAM_df.index = pd.MultiIndex.from_product([self.regions_io, self.sectors_io],
                                                       names=['region', 'sector'])

            # we consider that inputs of food in exiobase sectors are due to canteens and such -> exclude as we consider
            # employees would eat another way if canteen was not there -> canteen not required for the manufacture for
            # the commodity/service
            remove_canteen = pd.read_excel(pkg_resources.resource_filename(
                __name__, '/Data/ecoinvent/ei'+self.ei_version+'/double_counting_correction/force_canteen_out.xlsx'),
                index_col=0)

            # converts ecoinvent names to their corresponding exiobase sector in a dataframe format
            eco_to_exio = pd.get_dummies(
                self.filter['Hybridized processes'].set_index('code').loc[:, 'exiobase_sector']).astype(int).T
            multiindex = pd.MultiIndex.from_product([self.regions_io, eco_to_exio.index])
            eco_to_exio = pd.concat([eco_to_exio] * len(self.regions_io))
            eco_to_exio.index = multiindex
            eco_to_exio = eco_to_exio.reindex(self.A_io.index).fillna(0)

            # gamma represents the inputs that are deemed missing on purpose based on STAM
            gamma_filter_matrix = STAM_df.dot((STAM_table.mul(remove_canteen))).dot(
                STAM_df.T.dot(eco_to_exio.reindex(self.A_io.index).fillna(0)))
            gamma_filter_matrix /= len(self.regions_io)

            # save RAM
            del STAM_df
            del eco_to_exio
            del STAM_table
            del remove_canteen
            del self.A_io

            # instead of phi, we do it directly into lambda, saves RAM
            # we correct lambda for the presence of inputs of similar functionality within the LCA description
            # presence of one of these inputs (e.g., diesel) forces the value zero to all other functionally similar
            # EEIO inputs (e.g., gasoline or kerosene)
            for category in ['Liquid Fuels', 'Solid Fuels', 'Gaseous Fuels', 'Electricity/heat', 'Transport']:
                ix = lambda_filter_matrix.loc(axis=0)[:, io_categories[category]].index
                cols_with_zero = (lambda_filter_matrix.loc[ix] == 0).any(axis=0)
                lambda_filter_matrix.loc[ix, cols_with_zero] = 0

            self.logger.info("Applying the various double counting correction matrices to the uncorrected upstream "
                             "cutoff matrix...")

            # apply various filters to remove the double counting from the uncorrected upstream cutoffs matrix
            self.A_io_f_uncorrected = gamma_filter_matrix.multiply(
                lambda_filter_matrix.multiply(self.A_io_f_uncorrected))

            # we did the modifications on the uncorrected matrix to save RAM, now we rename
            self.A_io_f = self.A_io_f_uncorrected
            # and delete to save RAM
            del self.A_io_f_uncorrected

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

    def import_hybridized_database(self, cutoff=0.00001):
        """Function imports the hybridized version of ecoinvent into brightway2."""

        # get codes of exiobase in dictionary for speed
        db_exiobase_name = 'exiobase'
        exio_find_code_easily = {(i.as_dict()['location'], i.as_dict()['name']): i.as_dict()['code'] for i in
                                 bw2.Database(db_exiobase_name)}

        if cutoff != 0:
            self.logger.info("Remove irrelevant upstream cutoffs inputs...")

            # store prices in a df to use them to select irrelevant inputs
            prices = self.filter['Hybridized processes'].set_index('code').loc[:, 'price']

            # we want to remove some of the inputs because it takes a long time to write the database otherwise
            # so if an input is smaller than cutoff of the total price of the product -> irrelevant input
            for process_code in tqdm(list(self.A_io_f.columns), leave=True):
                self.A_io_f.loc[:, process_code] = self.A_io_f.loc[:, process_code].where(
                    self.A_io_f.loc[:, process_code] > prices.loc[process_code] * cutoff).fillna(0)

        self.logger.info("Reformat the upstream cutoffs matrix...")
        # change to pd.Series and remove null inputs to speed things up
        self.A_io_f = self.A_io_f.stack(future_stack=True).stack(future_stack=True)
        self.A_io_f = self.A_io_f[self.A_io_f != 0]
        # change to dict format to speed things up
        self.A_io_f = self.A_io_f.to_dict()

        self.logger.info("Copying existing regioinvent database...")
        hybrid_ei_data = bw2.Database(self.db_name).load()

        self.logger.info("Add relevant upstream cutoffs inputs to ecoinvent processes...")
        # put the upstream cutoffs in the ecoinvent descriptions
        for exc in tqdm(list(self.A_io_f.keys()), leave=True):
            # get the key of the ecoinvent process to hybridize
            receiving_process = (exc[3], exc[2])
            # get the key of the upstream cutoff to add to the ecoinvent process
            sectoral_input_bought = exio_find_code_easily[(exc[0], exc[1])]
            # add the exchange to the dictionary
            hybrid_ei_data[receiving_process]['exchanges'].append({
                "amount": self.A_io_f[exc],
                "name": exc[1],
                "location": exc[0],
                "type": "technosphere",
                "unit": "euro",
                "flow": str(uuid.uuid4()),
                "input": (db_exiobase_name, sectoral_input_bought)
            })

            # Remove the processed key to free memory
            del self.A_io_f[exc]

        name_hybrid_db = 'hybrid ' + self.db_name

        self.logger.info("Change name of the ecoinvent database to hybrid ecoinvent database...")
        # change the keys of processes from "ecoinvent" to "hybrid ecoinvent"
        hybrid_ei_data = {(name_hybrid_db, k[1]): v for k, v in hybrid_ei_data.items()}
        # change the keys of inputs and outputs from "ecoinvent" to "hybrid ecoinvent"
        for k in tqdm(hybrid_ei_data.keys(), leave=True):
            hybrid_ei_data[k]['database'] = name_hybrid_db
            for exc in hybrid_ei_data[k]['exchanges']:
                if exc['type'] in ['technosphere', 'production']:
                    if exc['input'][0] == self.db_name:
                        exc['input'] = (name_hybrid_db, exc['input'][1])
                    # for exiobase database, no output key was explicitly defined
                    if 'output' in exc.keys():
                        if exc['output'][0] == self.db_name:
                            exc['output'] = (name_hybrid_db, exc['output'][1])

        self.logger.info("Writing the hybrid database to brightway...")
        db = bw2.Database(name_hybrid_db)
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


class Hybridize_regioinvent:
    def __init__(self, bw2_project_name, regioinvent_db_name, ecoinvent_db_name, ecoinvent_version, path_to_exiobase):
        """
        Class hybridizes regioinvent with exiobase.

        :param bw2_project_name: [str] the name of your brightway2 project where you have an ecoinvent database stored
                                    that you wish to hybridize
        :param regioinvent_db_name: [str] the name of the regioinvent database you wish to hybridize
        :param ecoinvent_db_name: [str] the name of the ecoinvent database that was used as the basis for regioinvent
        :param ecoinvent_version: [str] the version of your ecoinvent database used as the basis for regioinvent.
                                        Supported versions = ['3.9', '3.9.1', '3.10', '3.10.1']
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
        self.db_ei_name = ecoinvent_db_name
        # check regioinvent database exists in project
        if regioinvent_db_name not in bw2.databases:
            raise KeyError(
                "The regioinvent database name passed does not match with any existing database within your brightway2 "
                "project.")
        self.db_name = regioinvent_db_name
        if ecoinvent_version not in ['3.9', '3.9.1', '3.10', '3.10.1']:
            raise KeyError(
                "Your ecoinvent version is not supported by the current version of pyLCAIO. Supported versions are "
                "'3.9', '3.9.1', '3.10' and '3.10.1'")
        if ecoinvent_version in ['3.9', '3.9.1']:
            self.ei_version = '3.9'
        elif ecoinvent_version in ['3.10', '3.10.1']:
            self.ei_version = '3.10'

        # load exiobase through pymrio
        self.logger.info("Loading exiobase...")
        io = pymrio.parse_exiobase3(path_to_exiobase)
        # we do not store the whole io object as many of the attributes are useless to us and we are short on RAM
        self.x_io = io.x.copy()
        self.A_io = io.A.copy()
        self.sectors_io = io.get_sectors().tolist()
        self.regions_io = io.get_regions().tolist()
        self.io_units = io.satellite.unit.copy()
        self.io_units[self.io_units == 'M.EUR'] = 'euro'
        # get the reference year of used exiobase version
        self.reference_year_IO = int(io.meta.description[-4:])
        # save RAM
        del io

        # load capital matrix of exiobase (capital endogenization)
        with gzip.open(pkg_resources.resource_filename(
                __name__, '/Data/capitals_exiobase/K_cfc_pxp_exio3.8.2_'+str(self.reference_year_IO)+'.gz.pickle'),
                'rb') as file:
            K = pickle.load(file)
        # add capital matrix to technology matrix. No need to adjust satellite accounts because they're not used
        self.A_io += K

        # load file with all filter and concordance information
        self.filter = pd.read_excel(pkg_resources.resource_filename(
            __name__, '/Data/regioinvent/ei'+self.ei_version+'/mappings/filters.xlsx'), None)
        self.filter_ei = pd.read_excel(pkg_resources.resource_filename(
            __name__, '/Data/ecoinvent/ei'+self.ei_version+'/mappings/filters.xlsx'), None)

        # concordance between regioinvent countries/regions and exiobase countries/regions
        with open(pkg_resources.resource_filename(
                __name__, '/Data/regioinvent/ei'+self.ei_version+'/mappings/geo_conc.json'), 'r') as f:
            self.concordance_geos = json.load(f)

        self.covered_geos = None
        self.not_capital_goods_codes = None
        self.codes_to_names = None
        self.A_io_f = None
        self.A_io_f_uncorrected = None

        self.get_relevant_info()

    def get_relevant_info(self):
        """ wurst is used to fasten the identification of processes with brightway2"""

        self.logger.info("Formatting relevant data of regioinvent for quick access...")

        # identify which geographies are covered for which product (useful for determining which countries are in RoW)
        self.covered_geos = {act.as_dict()['reference product']: [] for act in bw2.Database(self.db_name)}

        for act in bw2.Database(self.db_name):
            if act.as_dict()['location'] not in self.covered_geos[act.as_dict()['reference product']]:
                self.covered_geos[act.as_dict()['reference product']].append(act.as_dict()['location'])

        # identify which products in ecoinvent are not capital goods, useful for removing irrelevant inputs later
        self.not_capital_goods_codes = [act.as_dict()['code'] for act in bw2.Database(self.db_ei_name) if
                                        act.as_dict()['unit'] != 'unit']

        # get dictionary to quickly get a product name from the process code
        self.codes_to_names = {act.as_dict()['code']: act.as_dict()['reference product'] for act in
                               bw2.Database(self.db_name)}

    def get_uncorrected_upstream_cutoff_matrix(self):
        """
        Obtain the uncorrected upstream cutoff matrix.

        Brief methodology description:
        1. Obtain the concordance between regioinvent processes and exiobase sectors and identify which processes of
        eregiinvent are hybridized (in the Excel files of the Data folder of pylcaio)
        2. Loop through the "to-hybridize" processes and apply concordances, both sectoral and geographical
        3. Through the concordance, copy the corresponding exiobase sectors and apply the price (corrected for inflation)
        to get the uncorrected upstream cutoffs
        """

        self.logger.info("Getting the uncorrected upstream cutoff matrix...")

        # build the concordance matrix (both sector and geo concordance)
        concordance_matrix = pd.DataFrame(0, index=pd.MultiIndex.from_product([self.regions_io, self.sectors_io]),
                              columns=self.filter['Hybridized processes'].code, dtype='float')

        # loop through the difference processes to-be-hybridized
        for col in tqdm(concordance_matrix.columns, leave=True):
            act = bw2.Database(self.db_name).get(col)
            # extract location and corresponding exiobase sector
            geo = act.as_dict()['location']
            sector = self.filter['Hybridized processes'].loc[
                self.filter['Hybridized processes'].code == act.as_dict()['code'], 'exiobase_sector'].iloc[0]
            # if regioinvent process location is in exiobase regions (US -> US)
            if geo in list(concordance_matrix.index.levels[0]):
                concordance_matrix.loc[(geo, sector), col] = 1
            # if it needs some mapping
            elif geo in self.concordance_geos.keys():
                concordance_matrix.loc[(self.concordance_geos[geo], sector), col] = 1
            # special case for the dynamic region of ecoinvent: RoW
            else:
                # apply the weighted average of total prod (x) for relevant countries in RoW region
                concordance_matrix.loc[:, col] = (self.x_io.loc(axis=0)[[i for i in concordance_matrix.index.levels[0] if
                                                             i not in [
                                                                 set(self.covered_geos[self.codes_to_names[col]]) - {
                                                                     'RoW'}]], sector] /
                                      self.x_io.loc(axis=0)[[i for i in concordance_matrix.index.levels[0] if
                                                             i not in [
                                                                 set(self.covered_geos[self.codes_to_names[col]]) - {
                                                                     'RoW'}]], sector].sum()).reindex(
                    concordance_matrix.index).loc[:, 'indout'].fillna(0)

        # TODO get dynamic inflation numbers
        inflation = 1.25

        # upstream cutoff is basically the product of concordance (both sectoral and geographical) and price
        self.A_io_f_uncorrected = self.A_io.dot(concordance_matrix *
                                                self.filter['Hybridized processes'].set_index('code').price *
                                                inflation)

        # add a multiindex level to columns, being the name of the ecoinvent database
        self.A_io_f_uncorrected.columns = pd.MultiIndex.from_arrays(
            [[self.db_name] * len(self.A_io_f_uncorrected.columns), self.A_io_f_uncorrected.columns])

        # save RAM
        del self.x_io

    def determine_covered_inputs(self):
        """
        To correct for double counting te first step is to determine which inputs are already covered by the LCA
        database. To do so, we need to propagate the inputs of non-hybridized processes into inventories of hybridized
        inventories. So that if the production of timber (hybridized) requires power sawing (non-hybridized) the power
        saw is ultimately associated in the production of timber inventory so that an input of Machinery is not added
        through hybridization.
        Also, once all inputs have been propagated to hybridized processes, we only keep relevant inputs. What do we
        mean by that? Well if you have a random 1e-11 kg input of steel somewhere in your LCA description is very unlikely
        to mean that steel was properly covered by LCA. So how do we determine those relevant inputs? Well we fixed them
        to 1e-7 to acccount for cases where really small quantities are actually relevant, i.e., the use of metals for
        coating purposes for example.
        For capital goods specifically, the threshold has been set lower to acccount for the small quantities in LCA
        descriptions due to the distribution over the long lifetime of capital goods. The thrshold is therefore 1e-16
        for capital goods.
        """

        self.logger.info("Determining the covered inputs for each process of the LCI database...")

        # we need an LCA object to get to the technosphere matrix of brightway2
        lca = bw2.LCA({bw2.Database(self.db_ei_name).random(): 1}, list(bw2.methods)[0])
        lca.lci()

        # go to dense matrix
        A_ff = pd.DataFrame(lca.technosphere_matrix.todense(), lca.activity_dict.keys(),
                            lca.activity_dict.keys())
        # save RAM
        del lca
        # change notation (from LCA to IO notations)
        A_ff = pd.DataFrame(np.eye(len(A_ff)), A_ff.index, A_ff.columns) - A_ff

        # only keep technosphere matrix with non-hybridized processes
        A_not_hyb = A_ff.copy('deep')
        A_not_hyb.loc(axis=0)[:, self.filter_ei['Hybridized processes'].code] = 0

        # only keep technosphere matrix with hybridized processes
        A_hyb = A_ff.copy('deep')
        not_hyb = (self.filter_ei['Market processes'].code.tolist() +
                   self.filter_ei['Internal and activities'].code.tolist() +
                   self.filter_ei['Empty and aggregated processes'].code.tolist())
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
        self.covered_inputs = A_hyb.dot(inverse)

        # need to save RAM
        del A_hyb
        del inverse

        # remove irrelevant inputs for non capital goods
        self.covered_inputs.loc(axis=0)[:, self.not_capital_goods_codes] = (
            self.covered_inputs.loc(axis=0)[:, self.not_capital_goods_codes].mask(
                self.covered_inputs.loc(axis=0)[:, self.not_capital_goods_codes] < 1e-7, 0))

        # remove irrelevant inputs for capital goods
        self.covered_inputs[self.covered_inputs < 1e-16] = 0

        # remove first level multiindex
        self.covered_inputs = self.covered_inputs.droplevel(0)

        # only keep non-null ecoinvent inputs (propagated ones became null)
        self.covered_inputs = self.covered_inputs.loc[self.covered_inputs.sum(1)[self.covered_inputs.sum(1) != 0].index]

        # convert ecoinvent inputs into exiobase sectoral inputs
        self.covered_inputs.index = [
            self.filter_ei['Hybridized processes'].set_index('code').loc[i, 'exiobase_sector'] for i in
            self.covered_inputs.index]

        # groupby the sectoral inputs
        self.covered_inputs = self.covered_inputs.groupby(self.covered_inputs.index).sum()

        # only keep hybridized processes as columns
        self.covered_inputs = self.covered_inputs.loc(axis=1)[:, self.filter_ei['Hybridized processes'].code]

        # convert codes to names
        self.covered_inputs.columns = [bw2.Database(i[0]).get(i[1]).as_dict()['reference product'] for i in
                                       self.covered_inputs.columns]

        # groupby columns
        self.covered_inputs = self.covered_inputs.groupby(axis=1, level=0).sum()

    def correct_double_counting(self, method='STAM'):
        """
        Correcting the occurring double counting.

        :param method: [str] The method for double counting correction can be either "binary" or "STAM".
                            Default option = "STAM"
        """

        self.logger.info("Getting the lambda matrix for double counting correction...")

        # change size of covered_inputs to match exiobase sectors
        multi_index = pd.MultiIndex.from_product([self.regions_io, self.covered_inputs.index])
        self.covered_inputs = pd.concat([self.covered_inputs] * len(self.regions_io), axis=0)
        self.covered_inputs.index = multi_index
        self.covered_inputs = self.covered_inputs.reindex(self.A_io.index).fillna(0)

        # only keep inputs that are NOT covered by ecoinvent (i.e., equal to zero)
        self.covered_inputs = self.covered_inputs.mask(self.covered_inputs > 0)
        # changes the 0s to 1s
        self.covered_inputs[self.covered_inputs == 0] = 1
        # NaNs are the places where covered inputs were, fill with 0s
        self.covered_inputs = self.covered_inputs.fillna(0)

        # get the lambda filter matrix for regioinvent
        lambda_filter_matrix = pd.DataFrame(0, self.covered_inputs.index, self.A_io_f_uncorrected.columns, dtype=float)
        hybridized = self.filter['Hybridized processes'].code.tolist()
        for col in tqdm(lambda_filter_matrix.columns, leave=True):
            if col[1] in hybridized:
                lambda_filter_matrix.loc[:, col] = self.covered_inputs.loc[:, self.codes_to_names[col[1]]]

        # save RAM
        del self.covered_inputs
        del hybridized

        if method == 'binary':
            self.logger.info("Applying the various double counting correction matrices to the uncorrected upstream "
                             "cutoff matrix...")
            self.A_io_f = self.A_io_f_uncorrected.multiply(lambda_filter_matrix)

        elif method == 'STAM':

            self.logger.info("Getting the gamma matrix for double counting correction...")

            # need to check if a missing input is legit or due to an omission
            STAM_table = pd.read_excel(pkg_resources.resource_filename(
                __name__, '/Data/regioinvent/ei'+self.ei_version+'/double_counting_correction/STAM_table.xlsx'),
                index_col=0)
            with open(pkg_resources.resource_filename(
                    __name__, '/Data/regioinvent/ei'+self.ei_version+'/double_counting_correction/STAM_categories.txt'),
                    'r') as file:
                io_categories = eval(file.read())

            # STAM_df converts STAM_table to a format with exiobase sectors as indexes
            STAM_df = pd.DataFrame(0, index=self.sectors_io, columns=STAM_table.columns)
            for col in STAM_df.columns:
                STAM_df.loc[io_categories[col], col] = 1
            # extend from 200 sectors to the 9800 sectors of exiobase
            STAM_df = pd.concat([STAM_df] * len(self.regions_io), axis=0)
            # add first level of multiindex (regions)
            STAM_df.index = pd.MultiIndex.from_product([self.regions_io, self.sectors_io],
                                                    names=['region', 'sector'])
            # we consider that inputs of food in exiobase sectors are due to canteens and such -> exclude as we consider
            # employees would eat another way if canteen was not there -> canteen not required for the manufacture for
            # the commodity/service
            remove_canteen = pd.read_excel(pkg_resources.resource_filename(
                __name__, '/Data/regioinvent/ei'+self.ei_version+'/double_counting_correction/force_canteen_out.xlsx'),
                index_col=0)

            # converts ecoinvent names to their corresponding exiobase sector in a dataframe format
            eco_to_exio = pd.get_dummies(
                self.filter['Hybridized processes'].set_index('code').loc[:, 'exiobase_sector']).astype(int).T
            multiindex = pd.MultiIndex.from_product([self.regions_io, eco_to_exio.index])
            eco_to_exio = pd.concat([eco_to_exio] * len(self.regions_io))
            eco_to_exio.index = multiindex
            eco_to_exio = eco_to_exio.reindex(self.A_io.index).fillna(0)

            # gamma represents the inputs that are deemed missing on purpose based on STAM
            gamma_filter_matrix = STAM_df.dot((STAM_table.mul(remove_canteen))).dot(
                STAM_df.T.dot(eco_to_exio.reindex(self.A_io.index).fillna(0)))
            gamma_filter_matrix /= len(self.regions_io)

            # save RAM
            del STAM_df
            del eco_to_exio
            del STAM_table
            del remove_canteen
            del self.A_io

            self.logger.info("Getting the phi matrix for double counting correction...")
            # instead of phi, we do it directly into lambda, saves RAM
            # we correct lambda for the presence of inputs of similar functionality within the LCA description
            # presence of one of these inputs (e.g., diesel) forces the value zero to all other functionally similar
            # EEIO inputs (e.g., gasoline or kerosene)
            for category in ['Liquid Fuels', 'Solid Fuels', 'Gaseous Fuels', 'Electricity/heat', 'Transport']:
                ix = lambda_filter_matrix.loc(axis=0)[:, io_categories[category]].index
                cols_with_zero = (lambda_filter_matrix.loc[ix] == 0).any(axis=0)
                lambda_filter_matrix.loc[ix, cols_with_zero] = 0

            self.logger.info("Applying the various double counting correction matrices to the uncorrected upstream "
                             "cutoff matrix...")

            # apply various filters to remove the double counting from the uncorrected upstream cutoffs matrix
            self.A_io_f_uncorrected = gamma_filter_matrix.multiply(lambda_filter_matrix.multiply(self.A_io_f_uncorrected))

            # we did the modifications on the uncorrected matrix to save RAM, now we rename
            self.A_io_f = self.A_io_f_uncorrected
            # and delete to save RAM
            del self.A_io_f_uncorrected
