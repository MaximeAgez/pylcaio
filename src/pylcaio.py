"""pyLCAIO: handle and hybridize LCA and EEIO matrices

This module defines three classes (DatabaseLoader, LCAIO and Analysis) that loads an LCA and IO database and their
parameters, to then hybridize these databases and analyze the results.

Dependencies
-------------

- numpy
- pandas
- collections
- uuid
- os
- gzip
- pickle
- ast
- pkg_resources

Though not strictly speaking a dependency, this module relies on the
functionality of a pyMRIO object for reading in the IO tables.
"""

import numpy as np
import pandas as pd
import scipy.io
from collections import defaultdict
import os
import gzip
import pickle
import ast
import pkg_resources
import time

pd.set_option('mode.chained_assignment', None)
# pylint: disable-msg=C0103


class DatabaseLoader:
    """ Loads the LCA and IO databases to hybridize along with their parameters.

    Object instance variables:
    -------------------------

            - lca_database_name_and_version : name and version of the LCA database
            - io_database_name_and_version  : name and version of the IO database
            - LCA_database  : LCA database in a dictionnary of pandas dataframe
            - IO_database   : IO database in matrix formats

            - PRO_f     : metadata of the LCA database processes
            - A_ff      : technology matrix of the LCA database
            - A_io      : technology matrix of the IO database (commodity by commodity)
            - A_io_f    : upstream cut-off technology matrix
            - F_f       : environmental matrix of the LCA database
            - F_io      : environmental matrix of the IO database
            - F_io_f    : upstream cut-off environmental matrix
            - Y_io      : final demand associated to the base year of the IO database
            - y_f       : functional units for the processes of the LCA database
            - C_f       : characterization factors of the LCA database
            - C_io      : characterization factors of the IO database
            - STR_f     : stressors metadata of the LCA database
            - STR_IO    : stressors metadata of the IO database

            - io_categories : dictionary matching categories of products to product groups of the IO database

            - listcountry   : list regrouping the countries used in the IO database
            - listregions   : list regrouping the unique regions used in the LCA database
            - countries_per_regions   : dict matching unique regions of the LCA database to the corresponding
                                          countries of the IO database

            - replacements1 : dictionary regrouping replacements of geographies needing to occur separately
            - replacements2 : dictionary regrouping replacements of geographies needing to occur separately
            - replacements3 : dictionary regrouping replacements of geographies needing to occur separately

            - reference_year_IO         : the reference year of the IO database
            - number_of_countries_IO    : the number of countries (excluding residual regions) of the IO database
            - number_of_RoW_IO          : the number of residual regions of the IO database
            - number_of_products_IO     : the number of product groups of the IO database

            - list_to_hyb       : list regrouping the UUIDs of LCA processes to hybridize
            - list_not_to_hyb   : list regrouping the UUIDs of LCA processes not to hybridize
            - listmarket        : list regrouping the UUIDs of market processes in the LCA database
            - dummyprocesses    : list regrouping the UUIDs of dummy processes of the LCA database
            - listnottransacted : list regrouping the UUIDs of processes of the LCA database not being transacted
            - null_price        : list regrouping the UUIDs of processes whose product have a null price
            - listguillotine    : list regrouping the UUIDs of processes whose quality is deemed insufficient

            - listcreated       : list regrouping the UUIDs of the processes added to the LCA database

    Object methods:
    ---------------
        * ecoinvent-exiobase_loader()
        * add_process_to_ecoinvent()
        * filter_productName_ecoinvent()

    """

    def __init__(self, lca_database_processed, io_database_processed,
                 lca_database_name_and_version='ecoinvent3.5', io_database_name_and_version='exiobase3'):
        """ Define DatabaseLoader object
        Args:
        -----
            * lca_database_processed        : the LCA database as a dictionary of pandas dataframes
            * io_database_processed         : the IO database as pandas dataframes
            * path_to_io_database           : the path to the IO database (on your hardware)
            * lca_database_name_and_version : the name and version of the LCA database to be hybridized
            * io_database_name_and_version  : the name and version of the IO database to be hybridized
        """

        self.lca_database_name_and_version = lca_database_name_and_version
        self.io_database_name_and_version = io_database_name_and_version

        self.PRO_f = pd.DataFrame()
        self.A_ff = pd.DataFrame()
        self.A_io = pd.DataFrame()
        self.A_io_f = pd.DataFrame()
        self.F_f = pd.DataFrame()
        self.F_io = pd.DataFrame()
        self.F_io_f = pd.DataFrame()
        self.y_io = pd.DataFrame()
        self.y_f = pd.DataFrame()
        self.C_f = pd.DataFrame()
        self.C_io = pd.DataFrame()
        self.STR_f = pd.DataFrame()
        self.STR_io = pd.DataFrame()

        self.K_io = pd.DataFrame()
        self.K_io_f = pd.DataFrame()

        self.io_categories = defaultdict(list)
        self.categories_same_functionality = []

        self.listcountry = []
        self.listregions = []
        self.countries_per_regions = defaultdict(list)

        self.replacements1 = {}
        self.replacements2 = {}
        self.replacements3 = {}

        self.reference_year_IO = 0
        self.number_of_countries_IO = 0
        self.number_of_RoW_IO = 0
        self.number_of_products_IO = 0

        self.list_to_hyb = []
        self.list_not_to_hyb = []
        self.listmarket = []
        self.dummyprocesses = []
        self.listnottransacted = []
        self.null_price = []
        self.listguillotine = []
        self.list_uncovered_geographies = []

        self.listcreated = []

        self.LCA_database = lca_database_processed
        self.IO_database = io_database_processed

        del lca_database_processed
        del io_database_processed

        versions_of_ecoinvent = ['ecoinvent3.5', 'ecoinvent3.3']
        versions_of_exiobase = ['exiobase2', 'exiobase3']
        if self.lca_database_name_and_version not in versions_of_ecoinvent:
            print('The LCA database version you entered is not supported currently')
        if self.io_database_name_and_version not in versions_of_exiobase:
            print('The IO database version you entered is not supported currently')

    def combine_ecoinvent_exiobase(self, path_to_io_database, path_to_capitals):
        """ Loads every needed parameter to hybridize ecoinvent with exiobase as well as both databases
        Args
        ---
            * path_to_io_database   : the path leading to the io database folder
        """

        version_ecoinvent = extract_version_from_name(self.lca_database_name_and_version)
        version_exiobase = extract_version_from_name(self.io_database_name_and_version)

        self.PRO_f = self.LCA_database['PRO'].copy()
        self.PRO_f.price = self.PRO_f.price.fillna(0)
        self.A_ff = self.LCA_database['A'].copy()
        self.A_ff = self.A_ff.astype(dtype='float32')
        self.A_io = self.IO_database.A.copy()
        self.A_io = self.A_io.astype(dtype='float32')
        self.A_io_f = pd.DataFrame(0, index=self.A_io.index, columns=self.A_ff.columns, dtype='float32')
        self.F_f = self.LCA_database['F'].copy()
        self.F_f = self.F_f.astype(dtype='float32')
        self.y_io = self.IO_database.Y.copy()
        self.y_io = self.y_io.astype(dtype='float32')
        self.C_f = self.LCA_database['C'].copy()
        self.C_f = self.C_f.astype(dtype='float32')
        self.STR_f = self.LCA_database['STR'].copy().drop('cas', axis=1)
        self.STR_f.columns = ['MATRIXID', 'FULLNAME', 'UNIT', 'comp', 'subcomp']

        self.number_of_products_IO = len([i for i in self.IO_database.get_sectors()])
        self.number_of_RoW_IO = 5
        self.number_of_countries_IO = len([i for i in self.IO_database.get_regions()]) - self.number_of_RoW_IO

        if version_exiobase == str(2):
            self.reference_year_IO = 2007
            self.IO_database.emissions.S.index = self.IO_database.emissions.S.index.tolist()
            self.IO_database.emissions.S.columns = self.IO_database.emissions.S.columns.tolist()
            self.F_io = pd.concat(
                [self.IO_database.emissions.S, self.IO_database.resources.S, self.IO_database.materials.S])
            self.F_io = self.F_io / 1000000
            for_update = self.IO_database.factor_inputs.S.loc[self.IO_database.factor_inputs.unit[
                self.IO_database.factor_inputs.unit != 'M.EUR'].dropna().index] / 1000000
            self.IO_database.factor_inputs.S.update(for_update)
            self.F_io = pd.concat([self.F_io, self.IO_database.factor_inputs.S])
            self.F_io = self.F_io.select_dtypes(include=['float']).apply(pd.to_numeric, downcast='float')
            c_emissions = pd.read_excel(path_to_io_database + 'characterisation_CREEA_version2.2.2.xlsx',
                                        'Q_emission')
            c_emissions.columns = self.IO_database.emissions.S.index
            c_emissions = c_emissions.drop(c_emissions.index[0])
            c_factorinputs = pd.read_excel(path_to_io_database + 'characterisation_CREEA_version2.2.2.xlsx',
                                           'Q_factorinputs')
            c_factorinputs.columns = self.IO_database.factor_inputs.S.index
            c_factorinputs = c_factorinputs.drop(c_factorinputs.index[0])
            c_materials = pd.read_excel(path_to_io_database + 'characterisation_CREEA_version2.2.2.xlsx',
                                        'Q_materials')
            c_materials.columns = self.IO_database.materials.S.index
            c_materials = c_materials.drop(c_materials.index[0])
            c_resources = pd.read_excel(path_to_io_database + 'characterisation_CREEA_version2.2.2.xlsx',
                                        'Q_resources')
            c_resources.columns = self.IO_database.resources.S.index
            c_resources = c_resources.drop(c_resources.index[0])
            c_emissions.index = c_emissions.index.values
            c_emissions.columns = c_emissions.columns.values
            c_resources.columns = c_resources.columns.values
            self.C_io = pd.concat([c_emissions, c_materials, c_resources, c_factorinputs], sort=False)
            self.C_io = self.C_io.fillna(0)
            self.C_io = self.C_io.astype(dtype='float32')
            self.STR_io = pd.DataFrame([self.IO_database.emissions.S.index.tolist(),
                                        [i[0] for i in self.IO_database.emissions.S.index],
                                        [i[1] for i in self.IO_database.emissions.S.index]],
                                       index=['MATRIXID', 'FULLNAME', 'comp']).transpose()
            self.STR_io.index = self.STR_io.MATRIXID
            self.STR_io = self.STR_io.drop('MATRIXID', axis=1)

        if version_exiobase == str(3):
            # removing digits in the product group names of exiobase 3
            new_index_tuple = []
            for index in self.IO_database.A.index:
                if any(char.isdigit() for char in index[1]):
                    new_index_tuple.append((index[0], index[1][:-5]))
                else:
                    new_index_tuple.append(index)
            index_without_numbers = pd.MultiIndex.from_tuples(new_index_tuple, names=['region', 'sector'])
            self.A_io.index = index_without_numbers
            self.A_io.columns = index_without_numbers
            self.y_io.index = index_without_numbers
            self.IO_database.calc_all()
            self.F_io = self.IO_database.satellite.S
            # emissions for millions of euros, we want them in euros, except value added
            for_update = self.IO_database.satellite.S.loc[self.IO_database.satellite.unit[
                self.IO_database.satellite.unit != 'M.EUR'].dropna().index]/1000000
            self.F_io.update(for_update)
            self.F_io.columns = index_without_numbers
            self.F_io = self.F_io.astype(dtype='float32')
            self.C_io = pd.concat([pd.read_excel(
                pkg_resources.resource_stream(__name__, '/Data/characterisation_CREEA_version3.xlsx'),
                'Q_emission'),
                pd.read_excel(
                    pkg_resources.resource_stream(__name__, '/Data/characterisation_CREEA_version3.xlsx'),
                    'Q_materials'),
                pd.read_excel(
                    pkg_resources.resource_stream(__name__, '/Data/characterisation_CREEA_version3.xlsx'),
                    'Q_resources'),
                pd.read_excel(
                    pkg_resources.resource_stream(__name__, '/Data/characterisation_CREEA_version3.xlsx'),
                    'Q_factor_inputs')], sort=False).fillna(0)
            self.C_io = self.C_io.astype(dtype='float32')
            self.reference_year_IO = int(self.IO_database.meta.description[-4:])

        self.F_io_f = pd.DataFrame(0, self.F_io.index, self.F_f.columns, dtype='float32')

        # CAPITAL GOODS
        K_dict = scipy.io.loadmat(path_to_capitals)
        Kbar = pd.DataFrame(K_dict['KbarCfc'].toarray())
        inv_diag_x = pd.DataFrame(np.diag(self.IO_database.x.iloc[:, 0]))
        for position in inv_diag_x:
            if inv_diag_x.loc[position, position] != 0:
                inv_diag_x.loc[position, position] = 1 / inv_diag_x.loc[position, position]

        self.K_io = Kbar.dot(inv_diag_x)
        self.K_io.index = self.A_io.index
        self.K_io.columns = self.A_io.columns
        self.K_io = self.K_io.astype('float32')

        del self.LCA_database
        del self.IO_database

        # STAM CATEGORIES

        self.io_categories = ast.literal_eval(pkg_resources.resource_string(__name__,
                '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
                                                                             '/STAM_categories.txt').decode(
            'utf-8'))
        self.categories_same_functionality = ast.literal_eval(
            pkg_resources.resource_string(
                __name__,
                '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/STAM_functional_categories.txt').decode(
                'utf-8'))

        # GEOGRAPHY CONCORDANCE

        self.listcountry = ast.literal_eval(
            pkg_resources.resource_string(
                __name__,
                '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/countries.txt').decode(
                'utf-8'))
        self.listregions = ast.literal_eval(
            pkg_resources.resource_string(
                __name__,
                '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/regions.txt').decode(
                'utf-8'))
        self.countries_per_regions = ast.literal_eval(
            pkg_resources.resource_string(
                __name__,
                '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
                '/countries_per_regions.txt').decode('utf-8'))
        self.replacements1 = ast.literal_eval(pkg_resources.resource_string(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
            '/geography_replacements_regions.txt').decode('utf-8'))
        self.replacements2 = ast.literal_eval(pkg_resources.resource_string(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
            '/geography_replacements_other.txt').decode('utf-8'))
        self.replacements3 = ast.literal_eval(pkg_resources.resource_string(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
            '/geography_replacements_RoW.txt').decode('utf-8'))

        self.PRO_f['io_geography'] = self.PRO_f.geography.copy()
        self.PRO_f.io_geography = self.PRO_f.io_geography.replace(self.replacements1, regex=True)
        self.PRO_f.io_geography = self.PRO_f.io_geography.replace(self.replacements2, regex=True)
        self.PRO_f.io_geography = self.PRO_f.io_geography.replace(self.replacements3, regex=True)
        # cannot replace these normally because they alter the existing regions of ecoinvent3.5
        if version_ecoinvent == str(3.5):
            self.PRO_f.io_geography[[
                i for i in self.PRO_f.index if self.PRO_f.io_geography[i] in ['ER', 'NA', 'TN']]] = 'WF'
            self.PRO_f.io_geography[
                [i for i in self.PRO_f.index if self.PRO_f.io_geography[i] in ['NI', 'AR']]] = 'WL'
            if version_exiobase == 2:
                self.PRO_f.io_geography[self.PRO_f.io_geography == 'HR'] = 'WE'

        # PRODUCT CONCORDANCE

        exceptions = pd.read_excel(pkg_resources.resource_stream(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Product_Concordances.xlsx'),
            'Exceptions')
        concordance = pd.read_excel(pkg_resources.resource_stream(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Product_Concordances.xlsx'),
            'Concordance per product')
        convert_sector_code = pd.read_excel(pkg_resources.resource_stream(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Product_Concordances.xlsx'),
            'Description_Exiobase')

        exceptions.index = exceptions.UUID
        concordance = concordance.drop(['activityName', 'productName'], 1)
        self.PRO_f = self.PRO_f.merge(concordance, 'outer')
        self.PRO_f.index = self.PRO_f.activityId + '_' + self.PRO_f.productId
        self.PRO_f.Concordance.update(exceptions.Concordance)
        # convert exiobase codes (e.g. p01) to names of the sectors (e.g. 'Paddy rice')
        self.PRO_f = self.PRO_f.merge(convert_sector_code, left_on='Concordance', right_on='EXIO_code', how='left')
        self.PRO_f = self.PRO_f.drop(['Concordance', 'EXIO_code'], axis=1)
        self.PRO_f.index = self.PRO_f.activityId + '_' + self.PRO_f.productId
        # We want the indexes of PRO_f and A_ff in the same order
        self.PRO_f = self.PRO_f.reindex(self.A_ff.index)

        # LOADING THE FILTER

        self.list_to_hyb = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Hybridized').index.tolist()
        self.listmarket = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Market').index.tolist()
        self.listnottransacted = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Not commercialized').index.tolist()
        self.listguillotine = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Poor quality').index.tolist()
        self.dummyprocesses = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Empty processes').index.tolist()
        self.null_price = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'No price').index.tolist()
        self.list_uncovered_geographies = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Uncovered geography').index.tolist()
        self.list_not_to_hyb = (
                self.listmarket + self.listnottransacted + self.listguillotine + self.dummyprocesses
                + self.null_price + self.list_uncovered_geographies)

        self.read_template()

        self.qualitychecks()

        return LCAIO(PRO_f=self.PRO_f, A_ff=self.A_ff, A_io=self.A_io, A_io_f=self.A_io_f, F_f=self.F_f, F_io=self.F_io,
                     F_io_f=self.F_io_f, y_io=self.y_io, C_f=self.C_f, C_io=self.C_io, STR_f=self.STR_f,
                     STR_io=self.STR_io, listcountry=self.listcountry, listregions=self.listregions, K_io=self.K_io,
                     countries_per_regions=self.countries_per_regions, reference_year_IO=self.reference_year_IO,
                     number_of_countries_IO=self.number_of_countries_IO, number_of_RoW_IO=self.number_of_RoW_IO,
                     number_of_products_IO=self.number_of_products_IO, list_to_hyb=self.list_to_hyb,
                     list_not_to_hyb=self.list_not_to_hyb, listmarket=self.listmarket,
                     dummyprocesses=self.dummyprocesses, listnottransacted=self.listnottransacted,
                     null_price=self.null_price, listguillotine=self.listguillotine,
                     list_uncovered_geographies=self.list_uncovered_geographies, io_categories=self.io_categories,
                     categories_same_functionality=self.categories_same_functionality,
                     lca_database_name_and_version=self.lca_database_name_and_version,
                     io_database_name_and_version=self.io_database_name_and_version)

    def read_template(self):
        template_foreground_metadata = template_sheet_treatment(pd.read_excel(pkg_resources.resource_stream(__name__,
                                                      '/Template.xlsx'), 'Metadata_foreground'))
        template_foreground_exchanges = template_sheet_treatment(pd.read_excel(pkg_resources.resource_stream(__name__,
                                                      '/Template.xlsx'), 'Unit_processes_exchanges')).ffill()

        if template_foreground_metadata.isna().all().all():
            return
        if template_foreground_exchanges.isna().all().all():
            return

        for new_process_to_hybridize in [i for i in template_foreground_metadata.index
                                         if template_foreground_metadata.to_hybridize[i] == 'yes']:
            self.list_to_hyb.append(new_process_to_hybridize)
        for new_process_not_to_hybridize in [i for i in template_foreground_metadata.index
                                         if template_foreground_metadata.to_hybridize[i] == 'no']:
            self.list_not_to_hyb.append(new_process_not_to_hybridize)
        for new_process_not_to_hybridize in [i for i in template_foreground_metadata.index
                                         if str(template_foreground_metadata.to_hybridize[i]) == 'nan']:
            self.list_not_to_hyb.append(new_process_not_to_hybridize)
        template_foreground_metadata.drop('to_hybridize', axis=1, inplace=True)
        self.PRO_f = pd.concat([self.PRO_f, template_foreground_metadata], sort=False)

        for process_to_add in template_foreground_metadata.index:
            self.add_process_to_matrices(process_to_add)

        technosphere_exchanges = template_foreground_exchanges.loc[[i for i in template_foreground_exchanges.index if
                                                                    template_foreground_exchanges.input_output[
                                                                        i] in self.PRO_f.index]]
        pivot1 = pd.pivot_table(technosphere_exchanges.loc[:, ['ProcessId', 'input_output', 'value']],
                            values='value', index=['input_output'], columns=['ProcessId'], aggfunc=sum)
        del pivot1.index.name
        del pivot1.columns.name
        pivot1 = self.LCA_convention_to_IO(pivot1)
        self.A_ff.loc[pivot1.index, template_foreground_metadata.index] = pivot1

        biosphere_exchanges = template_foreground_exchanges.loc[[i for i in template_foreground_exchanges.index if
                                                                    template_foreground_exchanges.input_output[
                                                                        i] in self.STR_f.index]]

        pivot2 = pd.pivot_table(biosphere_exchanges.loc[:, ['ProcessId', 'input_output', 'value']],
                                values='value', index=['input_output'], columns=['ProcessId'], aggfunc=sum)
        del pivot2.index.name
        del pivot2.columns.name
        pivot2 = self.LCA_convention_to_IO(pivot2)
        self.F_f.loc[pivot2.index, template_foreground_metadata.index] = pivot2

    def add_process_to_matrices(self, index):
        """ Creates an empty row and empty column for processes to add """
        self.A_ff[index] = 0
        self.A_ff = pd.concat([self.A_ff, pd.DataFrame(0, columns=self.A_ff.columns, index=[index])], sort=False)
        self.F_f[index] = 0
        # self.y_f[index] = 0
        # self.y_f.loc[index, index] = 1

    def LCA_convention_to_IO(self, dataframe):
        """ Changes the convetion of an LCA technology matrix from LCA to IO """
        # no ones on the diagonal
        dataframe[dataframe == 1] = 0
        # only positive values
        dataframe = dataframe.abs()
        return dataframe

    def qualitychecks(self):
        """ Several quality checks to ensure that data entered by the user in the numerous excel and text files
        does not contain errors """
        if self.A_ff.empty:
            raise Exception('A_ff is empty!')
        if self.A_io.empty:
            raise Exception('A_io is empty!')
        if self.F_f.empty:
            raise Exception('F_f is empty!')
        if self.F_io.empty:
            raise Exception('F_io is empty!')
        if self.C_f.empty:
            raise Exception('C_f is empty!')
        if self.C_io.empty:
            raise Exception('C_io is empty!')
        if len(self.list_to_hyb) == 0:
            raise Exception('list_to_hyb is empty!')
        if len(self.list_not_to_hyb) == 0:
            raise Exception('list_not_to_hyb is empty!')
        if len(self.listmarket) == 0:
            raise Exception('listmarket is empty!')
        if len(self.listcountry) == 0:
            raise Exception('listcountry is empty!')
        if len(self.listregions) == 0:
            raise Exception('listregions is empty!')
        if len(self.countries_per_regions) == 0:
            raise Exception('no concordance between regions and countries entered')
        if len(self.io_categories) == 0:
            raise Exception('no concordance between categories and product groups entered')
        if self.reference_year_IO == 0:
            raise Exception('no reference year entered')
        if self.number_of_products_IO == 0:
            raise Exception('no number of products of IO entered')
        if self.number_of_countries_IO == 0:
            raise Exception('no number of countries of IO entered')
        if len([k for k in self.countries_per_regions.keys() if k not in self.listregions]) != 0:
            raise Exception('countries_per_regions countains regions that are not in listregions')
        if len([k for k in self.listregions if k not in self.countries_per_regions.keys()]) != 0:
            raise Exception('some regions of listregions are not translated into countries in countries_per_regions')
        for list_of_country in self.countries_per_regions.values():
            for country in list_of_country:
                if country not in self.listcountry:
                    raise Exception('some countries in countries_per_regions are not present in listcountry')
        if len([cat for cat in self.categories_same_functionality if cat not in self.io_categories.keys()]) != 0:
            raise Exception('All categories of same functionality must also be declared in io_categories')
        groups = []
        for list_of_groups in self.io_categories.values():
            for group in list_of_groups:
                if group in self.A_io.index.levels[1]:
                    groups.append(group)
        if len(set(groups)) != len(groups):
            raise Exception('Each product group can only belong to one category')
        if pd.DataFrame(self.PRO_f.ProductTypeName[self.list_to_hyb]).isnull().any()[0]:
            raise Exception('A process to hybridize is not matched to any product group')
        if len([geo for geo in self.PRO_f.io_geography if (geo not in self.A_io.index.levels[0]
                                                       and geo not in self.listregions
                                                       and geo != 'RoW'
                                                      )]) != 0:
            raise Exception('A geography supposed to be used to hybridize processes is not included in the IO database')
        if not len(set(self.PRO_f.index.tolist())) == len(self.PRO_f.index.tolist()):
            raise Exception('Duplicate in the indexes of PRO_f')
        if self.A_ff.isnull().any().any():
            raise Exception('NaN values in A_ff')
        # if not self.A_ff[self.A_ff < 0].isna().all().all():
        #     raise Exception('Negative values in A_ff')
        if not (self.A_ff.index == self.PRO_f.index).all():
            raise Exception('A_ff and PRO_f do not have the same index')
        if not (self.A_ff.index == self.A_ff.columns).all():
            raise Exception('A_ff indexes and columns are not identical')
        if not len(self.list_to_hyb)+len(self.list_not_to_hyb) == len(self.PRO_f):
            raise Exception('Some processes are neither in list_to_hyb nor list_not_to_hyb')
        if not (self.F_f.index.tolist().sort() == self.C_f.columns.tolist().sort()):
            raise Exception('Rows of F_f must match with columns of C_f')
        if not (self.F_io.index.tolist().sort() == self.C_io.columns.tolist().sort()):
            raise Exception('Rows of F_io must match with columns of C_io')
        if not (self.A_ff.index.tolist().sort() == self.F_f.columns.tolist().sort()):
            raise Exception('Rows of A_ff must match with columns of F_f')
        if not (self.A_io.index.tolist().sort() == self.F_io.columns.tolist().sort()):
            raise Exception('Rows of A_io must match with columns of F_io')
        if not (self.PRO_f.activityId+'_'+self.PRO_f.productId == self.PRO_f.index).all():
            raise Exception('Indexes must be a combination of activityId and productId, in this order')


class LCAIO:
    """ Handles and hybridized LCA inventory matrices and EEIO tables

    Object instance variables and notation:
    --------------------------------------

        - H : Concordance matrix matching processes of the LCA database to product groups of the IO database
        - G : Concordance matrix matching product groups of the IO database to categories used in STAM

        - A_ff_processed        : the extended LCA database technology matrix where inputs of the processes not to
                                  hybridize are added to the unit process inventory of processes to hybridize
        - total_prod_country    : the volume production of each product group for each country of the IO database
        - total_prod_region     : the volume production of each product group for each region of the LCA database
        - total_prod_RoW        : the volume production of each product group for each RoW of the LCA database
        - dictRow               : the dictionary matching RoW regions to the countries they include
        - STAM_table            : the matrix based on heuristics used to enhance the correstion of double counting
        - A_io_f_uncorrected    : the non-corrected upstream cut-off matrix


        Key concatenated matrices:
            A     : all normalized product requirements (technical coefficients)
            F     : all normalized extensions
            C     : all characterisation factors
            y     : final demand of the hybrid system
            STR   : all stressors

    Object methods:
    ---------------
        * hybridize()
        * identify_rows()
        * calc_productions()
        * extend_inventory()
        * save_system()

    """

    def __init__(self, **kwargs):
        """ Define LCAIO object """

        self.lca_database_name_and_version = ''
        self.io_database_name_and_version = ''

        self.PRO_f = pd.DataFrame()
        self.A_ff = pd.DataFrame()
        self.A_io = pd.DataFrame()
        self.A_io_f = pd.DataFrame()
        self.F_f = pd.DataFrame()
        self.F_io = pd.DataFrame()
        self.F_io_f = pd.DataFrame()
        self.y_io = pd.DataFrame()
        self.C_f = pd.DataFrame()
        self.C_io = pd.DataFrame()
        self.STR_f = pd.DataFrame()
        self.STR_io = pd.DataFrame()

        self.K_io = pd.DataFrame()
        self.K_io_f = pd.DataFrame()

        self.io_categories = defaultdict(list)
        self.categories_same_functionality = []

        self.listcountry = []
        self.listregions = []
        self.countries_per_regions = defaultdict(list)

        self.replacements1 = {}
        self.replacements2 = {}
        self.replacements3 = {}

        self.reference_year_IO = 0
        self.number_of_countries_IO = 0
        self.number_of_RoW_IO = 0
        self.number_of_products_IO = 0

        self.list_to_hyb = []
        self.list_not_to_hyb = []
        self.listmarket = []
        self.dummyprocesses = []
        self.listnottransacted = []
        self.null_price = []
        self.listguillotine = []
        self.list_uncovered_geographies = []

        self.A_ff_processed = pd.DataFrame()
        self.total_prod_country = pd.DataFrame()
        self.total_prod_region = pd.DataFrame()
        self.total_prod_RoW = pd.DataFrame()
        self.dictRoW = {}
        self.STAM_table = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/STAM_table.xlsx')).fillna(0)
        self.H = pd.DataFrame()
        self.G = pd.DataFrame()
        self.A_io_f_uncorrected = pd.DataFrame()

        self.double_counting = ''
        self.capitals = ''
        self.description = []

        allowed_keys = list(self.__dict__.keys())
        self.__dict__.update((key, value) for key, value in kwargs.items() if key in allowed_keys)

    # ----------------------------CORE METHOD-------------------------------------

    def hybridize(self, method_double_counting, capitals_method=False):
        """
        Hybridize the LCA database with the IO database

        self.A_io_f_uncorrected is calculated following the equation (1) of the paper [insert doi]
        self.A_io_f is calculated following the euqation (2) of the paper [insert doi]

        Args:
        -----
            * method_double_counting: method to correct double counting with (='binary' or ='STAM')
        """

        if not method_double_counting:
            raise Exception('Please enter a method to correct double counting (i.e. binary or STAM)')

        if not capitals_method:
            self.description.append('Capitals not endogenized')
            print('Capitals will not be endogenized in the hybridization')

        self.identify_rows()
        self.update_prices_electricity()
        self.calc_productions()
        self.extend_inventory()

        self.A_io_f = pd.DataFrame(0.0, index=self.A_io.index, columns=self.A_ff.columns, dtype='float32')

        # ---- CONVERSION PART ------

        # product concordance matrix
        self.H = pd.DataFrame(0, index=self.A_io.index.get_level_values('sector').tolist()[
                                       0:self.number_of_products_IO], columns=self.A_ff.columns, dtype='int64')
        for sector in self.H.index:
            self.H.loc[sector, self.H.columns.intersection(self.PRO_f[self.PRO_f.ProductTypeName == sector].index)] = 1
        self.H = self.H.append([self.H] * (self.number_of_countries_IO + self.number_of_RoW_IO - 1))
        self.H.index = self.A_io.index

        # product concordance matrix filtered
        H_for_hyb = self.H.copy()
        H_for_hyb.loc[:, self.list_not_to_hyb] = 0

        # translate the geography concordance txt files into matrices
        matrix_countries_per_region = pd.DataFrame(0, index=self.listcountry,
                                                   columns=self.listcountry + self.listregions + list(
                                                       self.dictRoW.keys()))
        for country in self.listcountry:
            matrix_countries_per_region.loc[country, country] = 1
        for region in self.countries_per_regions.keys():
            matrix_countries_per_region.loc[self.countries_per_regions[region], region] = 1
        for RoW in self.dictRoW.keys():
            matrix_countries_per_region.loc[self.dictRoW[RoW], RoW] = 1

        # identify which region corresponds to which process
        region_covered_per_process = pd.DataFrame(0,
                                                  index=self.listcountry + self.listregions + list(self.dictRoW.keys()),
                                                  columns=self.A_ff.columns)
        for country in region_covered_per_process.index:
            df = pd.DataFrame(self.PRO_f.io_geography)[
                pd.DataFrame(self.PRO_f.io_geography) == country].dropna()
            region_covered_per_process.loc[country, region_covered_per_process.columns.intersection(df.index)] = 1

        # translate the previous region based coverage into country based coverage
        countries_covered_per_process = matrix_countries_per_region.dot(region_covered_per_process)
        countries_covered_per_process = pd.concat([countries_covered_per_process] * self.number_of_products_IO)
        countries_covered_per_process = countries_covered_per_process.sort_index()
        # annoying indexing processing
        listindex = []
        for country_position in range(0, len(countries_covered_per_process), self.number_of_products_IO):
            for commodity_position in range(0, self.number_of_products_IO):
                listindex.append((countries_covered_per_process.index[country_position],
                                  self.A_io.index.get_level_values('sector').tolist()[
                                  :self.number_of_products_IO][commodity_position]))
        concordance_geography = pd.concat([matrix_countries_per_region] * self.number_of_products_IO).sort_index()
        concordance_geography.index = listindex
        self.H.index = pd.Series(self.H.index.get_values())
        concordance_geography = concordance_geography.reindex(self.H.index)
        concordance_geography.index = self.A_io.index
        self.H.index = self.A_io.index

        # introduce production volumes in the mix
        concordance_geography_with_production_volumes = concordance_geography.multiply(pd.concat([self.total_prod_country] *
                                                                                  len(region_covered_per_process),
                                                                                  axis=1).values)
        # pivot the total_prod_country multiindex serie into a single index dataframe
        pivoting = pd.pivot_table(data=self.total_prod_country, values=self.total_prod_country,
                                  columns=self.total_prod_country.index.get_level_values('region'),
                                  index=self.total_prod_country.index.get_level_values('sector'))
        # some more indexing
        pivoting = pivoting.reindex(self.total_prod_country.index.get_level_values('sector')[:self.number_of_products_IO])
        pivoting.columns = pivoting.columns.droplevel(0)

        # add total prods of regions and RoWs
        production_volumes = pivoting.join(self.total_prod_region.join(self.total_prod_RoW))

        # previous dataframes were only for the number of products of IO and not for all the sector (=region+commodity)
        production_volumes_scaled_to_IO_sectors = pd.concat(
            [production_volumes] * (self.number_of_countries_IO + self.number_of_RoW_IO))
        production_volumes_scaled_to_IO_sectors.index = concordance_geography_with_production_volumes.index

        # weighted average with regards to production volumes
        weighted_concordance_geography = (
                    concordance_geography_with_production_volumes / production_volumes_scaled_to_IO_sectors).fillna(0)
        weighted_concordance_geography = weighted_concordance_geography.transpose().reindex(
            concordance_geography_with_production_volumes.columns).transpose()

        # inflation rate to consider the discrepancy between LCA database and IO database reference years
        inflation = get_inflation(self.reference_year_IO)

        Geo = weighted_concordance_geography.dot(region_covered_per_process)

        # the uncorrected A_io_f matrix is then defined by the following operation
        self.A_io_f_uncorrected = self.A_io.dot(H_for_hyb * inflation * Geo) * self.PRO_f.price
        self.A_io_f_uncorrected = self.A_io_f_uncorrected.astype(dtype='float32')

        # ------ DOUBLE COUNTING PART -------

        if method_double_counting == 'binary':
            lambda_filter_matrix = self.H.dot(self.A_ff_processed)
            lambda_filter_matrix = lambda_filter_matrix.mask(lambda_filter_matrix > 0)
            lambda_filter_matrix[lambda_filter_matrix == 0] = 1
            lambda_filter_matrix = lambda_filter_matrix.fillna(0)
            lambda_filter_matrix = lambda_filter_matrix.astype(dtype='float32')

            self.A_io_f = self.A_io_f_uncorrected.multiply(lambda_filter_matrix)

            self.double_counting = 'binary'
            self.description.append('binary')

        elif method_double_counting == 'STAM':
            lambda_filter_matrix = self.H.dot(self.A_ff_processed)
            lambda_filter_matrix = lambda_filter_matrix.mask(lambda_filter_matrix > 0)
            lambda_filter_matrix[lambda_filter_matrix == 0] = 1
            lambda_filter_matrix = lambda_filter_matrix.fillna(0)
            lambda_filter_matrix = lambda_filter_matrix.astype(dtype='float32')

            self.G = pd.DataFrame(0, index=self.A_io.index.get_level_values('sector').tolist()[
                                             :self.number_of_products_IO],
                                    columns=self.STAM_table.columns)
            for columns in self.G.columns:
                self.G.loc[[i for i in self.G.index if i in self.io_categories[columns]], columns] = 1
            self.G = self.G.append([self.G] * (self.number_of_countries_IO + self.number_of_RoW_IO - 1))
            self.G.index = self.A_io_f.index
            self.G.index = self.A_io_f.index

            gamma_filter_matrix = self.G.dot(self.STAM_table.dot(self.G.transpose().dot(self.H)))
            gamma_filter_matrix[gamma_filter_matrix == self.number_of_countries_IO + self.number_of_RoW_IO] = 1

            phi_filter_matrix = pd.DataFrame(1, index=self.A_io_f.index, columns=self.A_io_f.columns)
            categories_used_by_processes = self.G.transpose().dot(self.H.dot(self.A_ff_processed))
            categories_used_by_processes = categories_used_by_processes.astype(dtype='float32')

            for category in self.categories_same_functionality:
                list_category_to_zero = [i for i in categories_used_by_processes.columns if
                                         categories_used_by_processes.loc[category, i] != 0]
                phi_filter_matrix.loc[[i for i in phi_filter_matrix.index if
                                        i[1] in self.io_categories[category]], list_category_to_zero] = 0

            self.A_io_f = phi_filter_matrix.multiply(
                gamma_filter_matrix.multiply(lambda_filter_matrix.multiply(self.A_io_f_uncorrected)))
            self.A_io_f = self.A_io_f.astype(dtype='float32')

            self.double_counting = 'STAM'
            self.description.append('STAM')

        if capitals_method == 'LCA':
            self.description.append('Capitals endogenized prioritizing LCA capitals over IO')
            self.capitals = 'LCA'
            pass

        elif capitals_method == 'IO':
            # We use the concordance to determine which processes of ecoinvent are capitals
            capitals = ['Cattle', 'Construction work', 'Electrical machinery and apparatus n.e.c.',
                        'Machinery and equipment n.e.c.',
                        'Meat animals nec', 'Motor vehicles, trailers and semi-trailers',
                        'Office machinery and computers',
                        'Other transport equipment', 'Pigs', 'Poultry',
                        'Radio, television and communication equipment and apparatus',
                        'Medical, precision and optical instruments, watches and clocks',
                        'Furniture; other manufactured goods n.e.c.',
                        'Real estate services', 'Computer and related services']
            # K_io_f is defined the same as A_io_f
            self.K_io_f = self.K_io.dot(H_for_hyb * inflation * Geo) * self.PRO_f.price
            self.K_io_f = self.K_io_f.astype('float32')
            # the correction for double counting implies here to only trust IO capital and not LCA's
            self.A_ff.loc[[i for i in self.PRO_f.index if self.PRO_f.ProductTypeName[i] in capitals], :] = 0
            # Since we endogenized capitals we need to remove them from both final demand and factors of production
            self.y_io.loc[:, [i for i in self.y_io.columns if i[1] == 'Gross fixed capital formation']] = 0
            self.F_io.loc['Operating surplus: Consumption of fixed capital', :] = 0

            self.capitals = 'IO'
            self.description.append('Capitals endogenized prioritizing IO capitals over LCA')

        elif not capitals_method:
            self.capitals = False
            self.description.append('Capitals were not endogenized')

    # ---------------------PREPARATIONS FOR THE HYBRIDIZATION----------------------

    def identify_rows(self):

        unique_activities_using_row = list(set(
            self.PRO_f.activityNameId[[i for i in self.PRO_f.index if self.PRO_f.io_geography[i] == 'RoW']].tolist()))
        RoW_activities = defaultdict(list)
        tupl = [i for i in zip(self.PRO_f.activityNameId.loc[[i for i in self.PRO_f.index if self.PRO_f.activityNameId[
            i] in unique_activities_using_row]],
                               self.PRO_f.io_geography.loc[[i for i in self.PRO_f.index if self.PRO_f.activityNameId[
                                   i] in unique_activities_using_row]])]
        for activity, geography in tupl:
            RoW_activities[activity].append(geography)
        # remove RoW
        RoW_activities = {k: [v1 for v1 in v if v1 != 'RoW'] for k, v in RoW_activities.items()}
        # delete from RoW_activities processes which had only RoW as geography and are thus empty now
        for key in [i for i in list(RoW_activities.keys()) if RoW_activities[i] == []]:
            del RoW_activities[key]
        # put every element to the same level (elements that are lists are transformed to lists of lists)
        for values in list(RoW_activities.values()):
            for i in range(0, len(values)):
                if values[i] in self.countries_per_regions.keys():
                    values[i] = self.countries_per_regions[values[i]]
        # for elements that are lists of lists stemming from the replacement of ['RER'] by [['AT','BE',...]],
        # add all of the together in a single list
        for keys in RoW_activities.keys():
            for item in RoW_activities[keys]:
                if isinstance(item, list):
                    RoW_activities[keys] = sum_elements_list(RoW_activities[keys])
        # remove duplicates inside the elements
        for keys in list(RoW_activities.keys()):
            RoW_activities[keys] = list(set(RoW_activities[keys]))
        # need to sort to identify duplicates whose elements would be ordered differently and thus be treated as not
        # duplicated
        for keys in RoW_activities.keys():
            RoW_activities[keys].sort()
        # identify the combination of countries that are NOT inside the residual of each process
        dictactrow = {}
        residual_geo_IO_to_remove = ['WA', 'WE', 'WF', 'WL', 'WM']
        for keys in RoW_activities.keys():
            dictactrow[keys] = list(set(self.listcountry) - set(RoW_activities[keys]) - set(residual_geo_IO_to_remove))
        unique_RoWs = []
        for keys in dictactrow.keys():
            if dictactrow[keys] not in unique_RoWs:
                unique_RoWs.append(dictactrow[keys])
        # create name for the values of the different RoW
        listname = []
        for i in range(0, len(unique_RoWs)):
            listname.append('RoW' + '(' + str(i) + ')')
        # put all of that in dictRoW
        for i in range(0, len(unique_RoWs)):
            self.dictRoW[listname[i]] = unique_RoWs[i]
        try:
            # if RoWs are empty because processes from ecoinvent are too described
            del [[k for k in self.dictRoW.keys() if len(self.dictRoW[k]) == 0][0]]
        except IndexError:
            pass
        for keys in dictactrow:
            for keys2 in self.dictRoW:
                if dictactrow[keys] == self.dictRoW[keys2]:
                    dictactrow[keys] = keys2
        RoW_matrix = pd.DataFrame(list(dictactrow.values()), index=list(dictactrow.keys()), columns=['RoW_geography'])
        self.PRO_f = self.PRO_f.merge(RoW_matrix, left_on='activityNameId', right_on=RoW_matrix.index, how='outer')
        self.PRO_f.index = self.PRO_f.activityId + '_' + self.PRO_f.productId
        self.PRO_f = self.PRO_f.reindex(self.A_ff.index)
        self.PRO_f.io_geography.update(self.PRO_f.RoW_geography[self.PRO_f.io_geography == 'RoW'])
        self.PRO_f = self.PRO_f.drop('RoW_geography', axis=1)
        # might be some RoW or empty lists left in PRO_f
        self.PRO_f.io_geography[self.PRO_f.io_geography == 'RoW'] = 'GLO'
        self.PRO_f.io_geography.loc[[i for i in self.PRO_f.index if type(self.PRO_f.io_geography[i]) == list]] = 'GLO'

    def update_prices_electricity(self):
        """ Specially for ecoinvent and exiobase, but cannot be ran in the DatabaseLoader class as it requires
        identifying RoW regions to fully calculates all new prices """

        version_ecoinvent = extract_version_from_name(self.lca_database_name_and_version)
        version_exiobase = extract_version_from_name(self.io_database_name_and_version)

        electricity_price = pd.read_excel( pkg_resources.resource_stream(__name__,
                                      '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
                                      '/Regionalized_electricity_prices.xlsx'))

        electricity_processes = self.PRO_f.price[[i for i in self.PRO_f.index if self.PRO_f.price[i] == 0.0977]]

        if electricity_processes.empty:
            return

        merged = self.PRO_f.loc[electricity_processes.index.values,
                                ['price', 'io_geography']].merge(electricity_price,left_on='io_geography',
                                                                 right_on=electricity_price.index,how='left')
        merged.index = electricity_processes.index.values
        merged = merged.drop(['price', 'io_geography'], axis=1)
        self.PRO_f.price.update(merged.prices)

    def calc_productions(self):
        """ Calculates the different total productions for either countries, regions or RoWs."""

        # the user needs to determine the total demand before being able to calculate productions
        listdrop = []

        y = self.y_io.copy()
        y = y.groupby(level='region', axis=1).sum()
        i = np.eye(len(self.A_io))

        # solve the Leontief's system
        X = pd.DataFrame(np.linalg.solve(i - self.A_io, y), self.A_io.index, y.columns)

        # need to reindex the matrix calculated because python orders them alphabetically
        listindex = []
        for i in range(0, len(self.A_io.columns), self.number_of_products_IO):
            listindex.append(self.A_io.columns[i][0])
        X = X.reindex(listindex, axis=1)

        absent_countries = {}
        for i in range(0, len(list(self.countries_per_regions.values()))):
            absent_country = [item for item in self.listcountry if
                              item not in list(self.countries_per_regions.values())[i]]
            absent_countries[list(self.countries_per_regions.keys())[i]] = absent_country

        self.total_prod_country = pd.DataFrame(X.sum(axis=1), columns=['production'])

        listmatrixxx = []
        listlisteee = []
        listdfff = []
        for i in range(0, len(absent_countries)):
            listmatrixxx.append('matrixxx' + str(i))
            listlisteee.append('listeee' + str(i))
            listdfff.append('dfff' + str(i))
        listact = []
        for i in range(0, self.number_of_products_IO):
            listact.append(self.total_prod_country.index[i][1])
        for i in range(0, len(list(absent_countries.values()))):
            listadd = []
            listmatrixxx[i] = self.total_prod_country.drop(list(absent_countries.values())[i], axis=0, level=0)
            for k in range(0, self.number_of_products_IO):
                somme = 0
                for j in range(0, len(listmatrixxx[i]), self.number_of_products_IO):
                    somme += listmatrixxx[i].iloc[j + k, 0]
                listadd.append(somme)
            listlisteee[i] = listadd
            listdfff[i] = pd.DataFrame(listlisteee[i], listact, [list(absent_countries.keys())[i]])
            self.total_prod_region = self.total_prod_region.join(listdfff[i], how='outer')

        # next step we will consider the rest-of-the-World geographies, so the user has to run 'identify_RoWs' first
        if len(self.dictRoW) == 0:
            print('You need to run "identify_rows" before calculating the productions')
            return

        listmatrixxxx = []
        listlisteeee = []
        listdffff = []
        for k in range(0, len(list(self.dictRoW.keys()))):
            listmatrixxxx.append('matrixxxx' + str(k))
            listlisteeee.append('listeeee' + str(k))
            listdffff.append('dfff' + str(k))
            listdrop = []
            for i in range(0, len(self.dictRoW)):
                listadd = []
                for j in range(0, len(self.listcountry)):
                    if self.listcountry[j] not in list(self.dictRoW.values())[i]:
                        listadd.append(self.listcountry[j])
                listdrop.append(listadd)

        for i in range(0, len(list(self.dictRoW.keys()))):
            listadd = []
            listmatrixxxx[i] = self.total_prod_country.drop(listdrop[i], axis=0, level=0)
            for k in range(0, self.number_of_products_IO):
                somme = 0
                for j in range(0, len(listmatrixxxx[i]), self.number_of_products_IO):
                    somme += listmatrixxxx[i].iloc[j + k, 0]
                listadd.append(somme)
            listlisteeee[i] = listadd
            listdffff[i] = pd.DataFrame(listlisteeee[i], listact, [list(self.dictRoW.keys())[i]])
            self.total_prod_RoW = self.total_prod_RoW.join(listdffff[i], how='outer')

    def extend_inventory(self):
        """" Method creating a new technology matrix for the LCA database in which the inputs of processes not to
        hybridize are added to the unit process inventories of each LCA process to hybridize """
        # matrix of non-hybridized processes
        ANT = self.A_ff.copy()
        ANT.loc[self.list_to_hyb] = 0
        # matrix of production processes
        Amarket = self.A_ff.copy()
        Amarket.loc[self.listmarket] = 0

        Identity = pd.DataFrame(np.eye(len(ANT)), index=ANT.index, columns=ANT.columns)
        df = pd.DataFrame(np.linalg.inv(Identity - ANT), index=ANT.index, columns=ANT.columns)
        self.A_ff_processed = Amarket.dot(df)

        # the inversion brings a lot of noise which falsely contributes to double counting which is why it is removed
        self.A_ff_processed[self.A_ff_processed < 10 ** -16] = 0

    # -------------------------- EXPORT RESULTS -----------------------------------

    def save_system(self):
        """ Export the hybridized database to either dataframe via pickle """

        hybrid_system = {'PRO_f': self.PRO_f, 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                             'F_f': self.F_f, 'F_io': self.F_io, 'F_io_f': self.F_io_f,
                             'C_f': self.C_f, 'C_io': self.C_io, 'K_io': self.K_io, 'K_io_f': self.K_io_f,
                             'list_to_hyb': self.list_to_hyb, 'description': self.description}

        if self.capitals == 'IO' and self.double_counting == 'STAM':
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_IO_capitals_STAM/')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_IO_capitals_STAM/'))
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_IO_capitals_STAM/__init__.py')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_IO_capitals_STAM/__init__.py'))
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_IO_capitals_STAM/hybrid_system.pickle')),
                           'wb') as f:
                pickle.dump(hybrid_system, f)
        elif self.capitals == 'IO' and self.double_counting == 'binary':
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_IO_capitals_binary/')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_IO_capitals_binary/'))
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_IO_capitals_binary/__init__.py')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_IO_capitals_binary/__init__.py'))
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_IO_capitals_binary/hybrid_system.pickle')),
                           'wb') as f:
                pickle.dump(hybrid_system, f)
        elif self.capitals == 'LCA' and self.double_counting == 'STAM':
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_LCA_capitals_STAM/')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_LCA_capitals_STAM/'))
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_LCA_capitals_STAM/__init__.py')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_LCA_capitals_STAM/__init__.py'))
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_LCA_capitals_STAM/hybrid_system.pickle')),
                           'wb') as f:
                pickle.dump(hybrid_system, f)
        elif self.capitals == 'LCA' and self.double_counting == 'binary':
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_LCA_capitals_binary/')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_LCA_capitals_binary/'))
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_LCA_capitals_binary/__init__.py')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_LCA_capitals_binary/__init__.py'))
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_LCA_capitals_binary/hybrid_system.pickle')),
                           'wb') as f:
                pickle.dump(hybrid_system, f)
        elif not self.capitals and self.double_counting == 'STAM':
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_no_capitals_STAM/')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_no_capitals_STAM/'))
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_no_capitals_STAM/__init__.py')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_no_capitals_STAM/__init__.py'))
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_no_capitals_STAM/hybrid_system.pickle')),
                           'wb') as f:
                pickle.dump(hybrid_system, f)
        elif not self.capitals and self.double_counting == 'binary':
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_no_capitals_binary/')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_no_capitals_binary/'))
            if not os.path.exists(pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                            self.lca_database_name_and_version + '_' +
                                                                            self.io_database_name_and_version +
                                                                            '_no_capitals_binary/__init__.py')):
                os.mkdir(pkg_resources.resource_filename(__name__, '/Databases/' + self.lca_database_name_and_version +
                                                         '_' + self.io_database_name_and_version +
                                                         '_no_capitals_binary/__init__.py'))
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_no_capitals_binary/hybrid_system.pickle')),
                           'wb') as f:
                pickle.dump(hybrid_system, f)


class Analysis:
    """ Analyzes the results of the hybridization

    Object instance variables and notation:
    --------------------------------------
        - lca_database_name_and_version : name and version of the LCA database to analyze
        - io_database_name_and_version  : name and version of the IO database to analyze
        - d : matrix of environmental impacts

    Object methods:
    --------------
        * calc_lifecycle()
        * contribution_analysis_LCA_processes()
        * contribution_analysis_direct_upstream_cutoffs()
        * contribution_analysis_total()
        * look_into_A()

    """

    def __init__(self, lca_database_name_and_version, io_database_name_and_version, method_double_counting,
                 capitals_method=False):

        self.lca_database_name_and_version = lca_database_name_and_version
        self.io_database_name_and_version = io_database_name_and_version

        if capitals_method == 'IO' and method_double_counting == 'STAM':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_IO_capitals_STAM/hybrid_system.pickle')),
                           'rb') as f:
                self.hybrid_system = pd.read_pickle(f)
            self.capitals = 'IO'
            self.double_counting = 'STAM'
        elif capitals_method == 'IO' and method_double_counting == 'binary':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_IO_capitals_binary/hybrid_system.pickle')),
                           'rb') as f:
                self.hybrid_system = pd.read_pickle(f)
            self.capitals = 'IO'
            self.double_counting = 'binary'
        elif capitals_method == 'LCA' and method_double_counting == 'STAM':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_LCA_capitals_STAM/hybrid_system.pickle')),
                           'rb') as f:
                self.hybrid_system = pd.read_pickle(f)
            self.capitals = 'LCA'
            self.double_counting = 'STAM'
        elif capitals_method == 'LCA' and method_double_counting == 'binary':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_LCA_capitals_binary/hybrid_system.pickle')),
                           'rb') as f:
                self.hybrid_system = pd.read_pickle(f)
            self.capitals = 'LCA'
            self.double_counting = 'binary'
        elif not capitals_method and method_double_counting == 'STAM':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_no_capitals_STAM/hybrid_system.pickle')),
                           'rb') as f:
                self.hybrid_system = pd.read_pickle(f)
            self.capitals = 'no'
            self.double_counting = 'STAM'
        elif not capitals_method and method_double_counting == 'binary':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_no_capitals_binary/hybrid_system.pickle')),
                           'rb') as f:
                self.hybrid_system = pd.read_pickle(f)
            self.capitals = 'no'
            self.double_counting = 'binary'

        self.PRO_f = self.hybrid_system['PRO_f']
        self.A_ff = self.hybrid_system['A_ff']
        self.A_io = self.hybrid_system['A_io']
        self.A_io_f = self.hybrid_system['A_io_f']
        self.F_f = self.hybrid_system['F_f']
        self.F_io = self.hybrid_system['F_io']
        self.F_io_f = self.hybrid_system['F_io_f']
        self.C_f = self.hybrid_system['C_f']
        self.C_io = self.hybrid_system['C_io']
        self.K_io = self.hybrid_system['K_io']
        self.K_io_f = self.hybrid_system['K_io_f']
        self.list_to_hyb = self.hybrid_system['list_to_hyb']
        self.description = self.hybrid_system['description']

        del self.hybrid_system

        self.index_A = []
        for index in self.A_ff.index:
            self.index_A.append(index)
        for index in list(self.A_io.index.values):
            self.index_A.append(index)

        self.index_F = []
        for index in self.F_f.index:
            self.index_F.append(index)
        for index in list(self.F_io.index):
            self.index_F.append(index)

        self.y = np.zeros((len(self.A_ff)+len(self.A_io_f), len(self.A_ff)))
        np.fill_diagonal(self.y, 1, wrap=False)

        self.d = pd.DataFrame()

        self.L_A = np.array(None)
        self.L_K = np.array(None)
        self.L_AK = np.array(None)

        self.GWP100_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/GWP.txt').decode('utf-8'))
        self.Acidification_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/Acid.txt').decode('utf-8'))
        self.Eutrophication_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/Eutro.txt').decode('utf-8'))
        self.HTox_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/HTox.txt').decode('utf-8'))

    @property
    def F(self):
        """ Normalized extensions for whole system"""
        f = np.concatenate(
            [np.concatenate([self.F_f.values, np.zeros((len(self.F_f), len(self.A_io)))], axis=1),
             np.concatenate([self.F_io_f.values, self.F_io.values], axis=1)], axis=0)
        return f

    @property
    def C(self):
        """ Characterisation factors for whole system """
        try:
            c = pd.concat([self.C_f, self.C_io], sort=False).fillna(0)
        except TypeError:
            c = pd.concat([self.C_f, self.C_io]).fillna(0)
            c = c.reindex(self.F.index)
        return c

    def calc_lifecycle(self, capitals=True, stage='impacts'):
        """ Simply calculates lifecycle production, emissions, or impacts

        Args
        ----
            * stage:    either 'production', 'emissions', or 'impacts'
                        determines what is being calculated
            * capitals: boolean

        Returns
        -------
            * either lifecycle production (x), or emissions (e) or impact (d)
              vector

        """

        if capitals:
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/Leontief_inverses/L_AK.pickle')),
                           'rb') as f:
                self.L_AK = pd.read_pickle(f)
            x = self.L_AK.dot(self.y)

            if stage == 'production':
                return x

            e = self.F.dot(x)
            if stage == 'emissions':
                return e

            if stage == 'impacts':
                d = self.C.values.dot(e)
                self.d = pd.DataFrame(d, self.C.index, self.A_ff.columns)
                print('Calculations done! Results are contained in self.d')

        if not capitals:

            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/Leontief_inverses/L_A.pickle')),
                           'rb') as f:
                self.L_A = pd.read_pickle(f)

            x = self.L_A.dot(self.y)

            if stage == 'production':
                return x

            e = self.F.dot(x)
            if stage == 'emissions':
                return e

            if stage == 'impacts':
                d = self.C.values.dot(e)
                self.d = pd.DataFrame(d, self.C.index, self.A_ff.columns)
                print('Calculations done! Results are contained in self.d')

    def contribution_analysis_LCA_processes(self, UUID, impact_category):
        """ Contribution analysis for the LCA part of the hybrid LCI only
        Args:
        ----
            * UUID:             UUID of the process to examine
            * impact_category:  impact category to examine (GWP, ODP, Acidification, Eutrophication, HTox)
            """

        if not impact_category:
            print('Please enter an impact category: GWP100, Acidification, Eutrophication or HTox')
            return

        if impact_category == 'GWP100':
            generic_impact_category_name = self.GWP100_CML2001
        if impact_category == 'Acidification':
            generic_impact_category_name = self.Acidification_CML2001
        if impact_category == 'Eutrophication':
            generic_impact_category_name = self.Eutrophication_CML2001
        if impact_category == 'Human toxicity' or impact_category == 'HTox':
            generic_impact_category_name = self.HTox_CML2001

        for i in range(0, len(self.C_f)):
            try:
                if self.C_f.index[i] in generic_impact_category_name:
                    name_category = self.C_f.index[i]
            except NameError:
                print('For impact categories, enter GWP100, Acidification, Eutrophication or HTox')

        Identity = pd.DataFrame(np.eye(len(self.A_ff)), self.A_ff.index, self.A_ff.columns)
        X = pd.DataFrame(np.linalg.solve(Identity - self.A_ff.fillna(0), pd.DataFrame(np.diag(self.A_ff.loc[:, UUID]))),
                         index=self.A_ff.index, columns=self.A_ff.columns)
        matrix = pd.DataFrame(0, index=self.F_f.index, columns=self.F_f.columns)
        matrix.loc[:, UUID] = self.F_f.loc[:, UUID]
        D = pd.DataFrame(pd.DataFrame(self.C_f.fillna(0).dot(self.F_f.dot(X) + matrix)).loc[name_category, :],
                         index=self.F_f.columns, columns=[name_category])
        listproduct = [x for x in self.PRO_f.productName.loc[D.index]]
        listactivity = [x for x in self.PRO_f.activityName.loc[D.index]]
        D = D.join(pd.DataFrame(listproduct, index=D.index, columns=['productName']))
        D = D.join(pd.DataFrame(listactivity, index=D.index, columns=['activityName']))
        cols = D.columns.tolist()
        cols = cols[1:] + cols[:1]
        D = D[cols]
        return D.sort_values(D.columns[2], ascending=False)

    def contribution_analysis_direct_upstream_cutoffs(self, UUID, impact_category):
        """ Contribution analysis for the IO complements of the hybrid LCI only"
        Args:
        ----
            * UUID:             UUID of the process to examine
            * impact_category:  impact category to examine (GWP, ODP, Acidification, Eutrophication, HTox)
            """

        if not impact_category:
            print('Please enter an impact category: GWP100, Acidification, Eutrophication or HTox')
            return

        if impact_category == 'GWP100':
            generic_impact_category_name = self.GWP100_CML2001
        if impact_category == 'Acidification':
            generic_impact_category_name = self.Acidification_CML2001
        if impact_category == 'Eutrophication':
            generic_impact_category_name = self.Eutrophication_CML2001
        if impact_category == 'Human toxicity' or impact_category == 'HTox':
            generic_impact_category_name = self.HTox_CML2001

        for i in range(0, len(self.C_io)):
            try:
                if [self.C_io.index[i]] in generic_impact_category_name:
                    name_category = [self.C_io.index[i]]
            except NameError:
                print('For impact categories, enter GWP100, Acidification, Eutrophication or HTox')

        Y = pd.DataFrame(0.0, self.A_io_f.index, self.A_io.columns)
        Identity = pd.DataFrame(np.eye(len(self.A_io)), index=self.A_io.index, columns=self.A_io.columns)
        for i in range(0, len(Y)):
            Y.iloc[i, i] = self.A_io_f.iloc[i, self.A_io_f.columns.get_loc(UUID)]
        X = pd.DataFrame(np.linalg.solve(Identity - self.A_io, Y), index=self.A_io_f.index, columns=Y.columns)
        D = pd.DataFrame(self.C_io.dot(self.F_io.dot(X)).loc[name_category, :]).transpose()
        return D.sort_values(D.columns[0], ascending=False)

    def contribution_analysis(self, type_of_analysis, UUID, impact_category):

        name_impact_categories = self.check_impact_category(impact_category)

        if type_of_analysis == 'total':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_' + self.capitals + '_capitals_' +
                                                                      self.double_counting + '/L_AK.pickle')),
                           'wb') as f:
                self.L_AK = pd.read_pickle(f)

            X = self.L_AK.dot(np.diag(np.concatenate([self.A_ff.values[:, self.PRO_f.index.get_loc(UUID)],
                                                      self.A_io_f.values[:, self.PRO_f.index.get_loc(UUID)] +
                                                      self.K_io_f.values[:, self.PRO_f.index.get_loc(UUID)]], axis=0)))

        elif type_of_analysis == 'on_capitals_only':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_' + self.capitals + '_capitals_' +
                                                                      self.double_counting + '/L_K.pickle')),
                           'wb') as f:
                self.L_K = pd.read_pickle(f)

            X = self.L_K.dot(np.diag(np.concatenate([self.A_ff.values[:, self.PRO_f.index.get_loc(UUID)],
                                                     self.K_io_f.values[:, self.PRO_f.index.get_loc(UUID)]], axis=0)))

        elif type_of_analysis == 'without_capitals':
            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_' + self.capitals + '_capitals_' +
                                                                      self.double_counting + '/L_A.pickle')),
                           'wb') as f:
                self.L_A = pd.read_pickle(f)

            X = self.L_A.dot(np.diag(np.concatenate([self.A_ff.values[:, self.PRO_f.index.get_loc(UUID)],
                                                     self.A_io_f.values[:, self.PRO_f.index.get_loc(UUID)]], axis=0)))

        else:
            print('Enter the type_of_analysis desided')
            return

        results = self.get_results(UUID, X, name_impact_categories)
        return results

    def check_impact_category(self, impact_category):

        if not impact_category:
            print('Please enter an impact category: GWP100, Acidification, Eutrophication or HTox')
            return

        if impact_category == 'GWP100':
            generic_impact_category_name = self.GWP100_CML2001
        if impact_category == 'Acidification':
            generic_impact_category_name = self.Acidification_CML2001
        if impact_category == 'Eutrophication':
            generic_impact_category_name = self.Eutrophication_CML2001
        if impact_category == 'Human toxicity' or impact_category == 'HTox':
            generic_impact_category_name = self.HTox_CML2001

        for i in range(0, len(self.C_f)):
            try:
                if self.C_f.index[i] in generic_impact_category_name:
                    name_category_LCA = self.C_f.index[i]
            except NameError:
                print('For impact categories, enter GWP100, Acidification, Eutrophication or HTox')
        for i in range(0, len(self.C_io)):
            try:
                if [self.C_io.index[i]] in generic_impact_category_name:
                    name_category_IO = [self.C_io.index[i]]
            except NameError:
                print('For impact categories, enter GWP100, Acidification, Eutrophication or HTox')
        list_to_return = [name_category_LCA, name_category_IO]
        return list_to_return

    def get_results(self, UUID, X, name_categories):

        matrix = self.F.copy()
        matrix[:, self.PRO_f.index.get_loc(UUID) + 1:] = 0
        matrix[:, :self.PRO_f.index.get_loc(UUID) - 1] = 0
        d = self.C.values.dot(self.F.dot(X) + matrix)
        D = pd.DataFrame(d, index=self.C.index, columns=self.index_A)
        df = pd.DataFrame(D.loc[name_categories[0], :]).join(D.loc[name_categories[1], :].transpose())
        df['total'] = pd.DataFrame(
            df.iloc[:, df.columns.get_loc(name_categories[0])] + df.iloc[:, df.columns.get_loc(name_categories[1][0])])
        listproduct = [x for x in self.PRO_f.productName] + [float('nan')] * len(self.A_io)
        listactivity = [x for x in self.PRO_f.activityName] + [float('nan')] * len(self.A_io)
        df = df.join(pd.DataFrame(listproduct, index=df.index, columns=['productName']))
        df = df.join(pd.DataFrame(listactivity, index=df.index, columns=['activityName']))
        df = df.sort_values(by='total', ascending=False)
        return df

    def calc_Leontief_inverses(self, matrix):

        if matrix == 'L_A':
            a = np.concatenate(
                [np.concatenate([self.A_ff.values, np.zeros((len(self.A_ff), len(self.A_io)))], axis=1),
                 np.concatenate([self.A_io_f.values, self.A_io.values], axis=1)],
                axis=0)
            invert = np.linalg.inv(np.eye(len(self.A_ff) + len(self.A_io)) - a)
            invert = invert.astype('float32')

            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_' + self.capitals + '_capitals_' +
                                                                      self.double_counting + '/L_A.pickle')),
                           'wb') as f:
                pickle.dump(invert, f)

        if matrix == 'L_K':
            k = np.concatenate(
                [np.concatenate([self.A_ff.values, np.zeros((len(self.A_ff), len(self.A_io)))], axis=1),
                 np.concatenate([self.K_io_f.values, self.A_io.values + self.K_io.values], axis=1)],
                axis=0)
            invert = np.linalg.inv(np.eye(len(self.A_ff) + len(self.A_io)) - k)
            invert = invert.astype('float32')

            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_' + self.capitals + '_capitals_' +
                                                                      self.double_counting + '/L_K.pickle')),
                           'wb') as f:
                pickle.dump(invert, f)

        if matrix == 'L_AK':
            ak = np.concatenate(
                [np.concatenate([self.A_ff.values, np.zeros((len(self.A_ff), len(self.A_io)))], axis=1),
                 np.concatenate([self.A_io_f.values + self.K_io_f.values, self.A_io.values + self.K_io.values],
                                axis=1)],
                axis=0)
            invert = np.linalg.inv(np.eye(len(self.A_ff) + len(self.A_io)) - ak)
            invert = invert.astype('float32')

            with gzip.open((pkg_resources.resource_filename(__name__, '/Databases/' +
                                                                      self.lca_database_name_and_version + '_' +
                                                                      self.io_database_name_and_version +
                                                                      '_' + self.capitals + '_capitals_' +
                                                                      self.double_counting + '/L_AK.pickle')),
                           'wb') as f:
                pickle.dump(invert, f)


def extract_version_from_name(name_database):

    for i in range(0, len(name_database)):
        try:
            if type(int(name_database[i])) == int:
                return name_database[i:]
        except ValueError:
            pass


def get_inflation(reference_year):
    """ Returns the inflation rate between the year 2005 (base year for ecoinvent prices) and the reference year of
    the used IO database"""

    if reference_year == 1995:
        inflation = 0.83
    elif reference_year == 1996:
        inflation = 0.84
    elif reference_year == 1997:
        inflation = 0.86
    elif reference_year == 1998:
        inflation = 0.87
    elif reference_year == 1999:
        inflation = 0.88
    elif reference_year == 2000:
        inflation = 0.9
    elif reference_year == 2001:
        inflation = 0.92
    elif reference_year == 2002:
        inflation = 0.94
    elif reference_year == 2003:
        inflation = 0.96
    elif reference_year == 2004:
        inflation = 0.98
    elif reference_year == 2005:
        inflation = 1
    elif reference_year == 2006:
        inflation = 1.02
    elif reference_year == 2007:
        inflation = 1.04
    elif reference_year == 2008:
        inflation = 1.08
    elif reference_year == 2009:
        inflation = 1.08
    elif reference_year == 2010:
        inflation = 1.10
    elif reference_year == 2011:
        inflation = 1.13
    else:
        inflation = 1

    return inflation


def sum_elements_list(liste):
    concatenated_list = []
    for i in range(0,len(liste)):
        if isinstance(liste[i], list):
            concatenated_list += liste[i]
        else:
            concatenated_list += [liste[i]]
    return concatenated_list


def template_sheet_treatment(dataframe):
    """ Removes empty rows and potential typos created Unnamed columns in the template """
    dataframe = dataframe.dropna(how='all', axis=0)
    dataframe = dataframe.drop([j for j in [i for i in dataframe.columns] if 'Unnamed' in j], axis=1)
    return dataframe
