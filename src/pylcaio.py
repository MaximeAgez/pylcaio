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
import re
import pkg_resources
import warnings
from brightway2 import *
from bw2agg.scores import add_all_unit_score_exchanges_and_cfs
import hashlib

pd.set_option('mode.chained_assignment', None)

# pylint: disable-msg=C0103


class DatabaseLoader:
    """ Loads the LCA and IO databases to hybridize along with their parameters.

    Object instance variables:
    -------------------------

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

            - replacements : dictionary linking geographies from ecoinvent to exiobase

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

    Object methods:
    ---------------
        * combine_ecoinvent_exiobase()

    """

    def __init__(self, lca_database_processed, io_database_processed,
                 lca_database_name_and_version, io_database_name_and_version):
        """ Define DatabaseLoader object
        Args:
        -----
            * lca_database_processed        : the LCA database as a dictionary of pandas dataframes
            * io_database_processed         : the IO database as pandas dataframes
            * lca_database_name_and_version : the name and version of the LCA database to be hybridized
            * io_database_name_and_version  : the name and version of the IO database to be hybridized
        """

        warnings.filterwarnings("ignore", message="divide by zero encountered in divide")

        self.lca_database_name_and_version = lca_database_name_and_version
        self.io_database_name_and_version = io_database_name_and_version

        self.PRO_f = pd.DataFrame()
        self.A_ff = pd.DataFrame()
        self.A_io = pd.DataFrame()
        self.F_f = pd.DataFrame()
        self.F_io = pd.DataFrame()
        self.y_io = pd.DataFrame()
        self.IMP = pd.DataFrame()
        self.C_f = pd.DataFrame()
        self.C_io = pd.DataFrame()
        self.C_f_regio = pd.DataFrame()
        self.C_io_regio = pd.DataFrame()
        self.STR_f = pd.DataFrame()
        self.K_io = pd.DataFrame()
        self.X_io = pd.DataFrame()
        self.F_io_regio = pd.DataFrame()
        self.F_f_regio = pd.DataFrame()

        self.io_categories = defaultdict(list)
        self.categories_same_functionality = []
        self.processes_in_order = []
        self.sectors_of_IO = []
        self.regions_of_IO = []
        self.Y_categories = []
        self.flows_of_IO = []
        self.impact_methods_IO = []
        self.extended_flows_names = []
        self.extended_impact_names_CML = []
        self.extended_impact_names_IW_exio = []
        self.extended_impact_names_IW_eco = []
        self.regionalized_impact_names_exio = []
        self.regionalized_impact_names_eco = []
        self.regionalized_flow_names_exio = []
        self.regionalized_flow_names_eco = []

        self.listcountry = []
        self.countries_per_regions = defaultdict(list)

        self.replacements = {}

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

        self.description = []

        self.LCA_database = lca_database_processed
        self.IO_database = io_database_processed

        self.LCA_database['A'] = self.LCA_database['A'].fillna(0)
        self.LCA_database['F'] = self.LCA_database['F'].fillna(0)
        self.LCA_database['C'] = self.LCA_database['C'].fillna(0)

        del lca_database_processed
        del io_database_processed

        versions_of_ecoinvent = ['ecoinvent3.8', 'ecoinvent3.7.1', 'ecoinvent3.6', 'ecoinvent3.5']
        if self.lca_database_name_and_version not in versions_of_ecoinvent:
            raise ValueError('The ecoinvent version you entered is not supported currently. Supported versions are: '
                             'ecoinvent3.8, ecoinvent3.7.1, ecoinvent3.6 and ecoinvent3.5.')
        if self.io_database_name_and_version.split('exiobase')[1][0] != '3':
            raise ValueError('The exiobase version you entered is not supported currently. Pylcaio only supports '
                             'exiobase3')

    def combine_ecoinvent_exiobase(self, path_to_capitals='', complete_extensions=False,
                                   impact_world=False, regionalized=False):
        """ Loads every needed parameter to hybridize ecoinvent with exiobase as well as both databases
        Args:
        ----
            * path_to_capitals      : the path leading to the capitals new matrix (only required if wanting to
                                      endogenize capitals) only available for EXIOBASE3
            * complete_extensions   : boolean (True to complete exiobase3 extensions based on USEEIO extension)
                                      only available for EXIOBASE3
            * impact_world          : boolean (True to use Impact world + method, False to get classic methods)
                                      only available for EXIOBASE3 and ecoinvent3.5
            * regionalized          : boolean (True to use regionalized impacts) / only available for EXIOBASE3 and ecoinvent3.5
                                      available regionalized impacts are (if marked by *, only available for ecoinvent):
                                      - Terrestrial acidification
                                      - Freshwater eutrophication
                                      - Marine eutrophication
                                      - Freshwater acidification
                                      - Land transformation, biodiversity*
                                      - Land occupation, biodiversity
                                      - Thermally polluted water*
                                      - Water scarcity

        """

        # EXTRACTING DATA FROM ECOINVENT AND EXIOBASE

        version_ecoinvent = self.lca_database_name_and_version.split('ecoinvent')[1]
        version_exiobase = self.io_database_name_and_version.split('exiobase')[1][0]

        # extracting processes metadata from ecoinvent
        self.PRO_f = self.LCA_database['PRO'].copy()
        self.processes_in_order = self.PRO_f.index.tolist()
        del self.LCA_database['PRO']
        # extracting intermediate exchanges from ecoinvent
        self.PRO_f.price = self.PRO_f.price.fillna(0)
        self.LCA_database['A'].values[self.LCA_database['A'].values < 0] *= -1
        self.A_ff = scipy.sparse.csr_matrix(self.LCA_database['A'].values)
        del self.LCA_database['A']
        # extracting intermediate exchanges from exiobase
        self.A_io = scipy.sparse.csr_matrix(self.IO_database.A)
        self.sectors_of_IO = self.IO_database.get_sectors().tolist()
        self.regions_of_IO = self.IO_database.get_regions().tolist()
        self.Y_categories = self.IO_database.get_Y_categories().tolist()
        # extracting environmental exchanges from ecoinvent
        self.F_f = scipy.sparse.csr_matrix(self.LCA_database['F'].values)
        del self.LCA_database['F']
        # extracting final demand from exiobase
        self.y_io = scipy.sparse.csr_matrix(self.IO_database.Y)
        # extracting characterization factors from ecoinvent
        self.C_f = scipy.sparse.csr_matrix(self.LCA_database['C'].values)
        del self.LCA_database['C']
        # extracting methods and elementary flows metadata from ecoinvent
        self.IMP = self.LCA_database['IMP']
        self.STR_f = self.LCA_database['STR'].copy().drop('cas', axis=1)
        self.STR_f.columns = ['MATRIXID', 'FULLNAME', 'UNIT', 'comp', 'subcomp']

        # useful intel for later
        self.number_of_products_IO = len([i for i in self.IO_database.get_sectors()])
        self.number_of_RoW_IO = 5
        self.number_of_countries_IO = len([i for i in self.IO_database.get_regions()]) - self.number_of_RoW_IO

        # cleaning up exiobase sector names
        self.sectors_of_IO = [re.sub(r' \(\d\d\)', '', i) for i in self.sectors_of_IO]
        # extracting calculated matrices from exiobase (X, L, etc.)
        self.IO_database.calc_all()
        self.X_io = scipy.sparse.csr_matrix(self.IO_database.x.values)
        # exiobase is in million euros, so we divide by 1,000,000
        self.IO_database.satellite.S[9:] /= 1000000
        self.F_io = scipy.sparse.csr_matrix(self.IO_database.satellite.S)
        del self.IO_database.satellite.S
        del self.IO_database.satellite.F
        del self.IO_database.Z
        self.flows_of_IO = self.IO_database.satellite.get_index().tolist()
        # compiling a characterization matrix for exiobase
        self.C_io = pd.concat([pd.read_excel(
            pkg_resources.resource_stream(__name__, '/Data/characterisationEXIOBASE3_adaptedFromEXIOBASE2.xlsx'),
            'Q_emission'),
            pd.read_excel(
                pkg_resources.resource_stream(__name__,
                                              '/Data/characterisationEXIOBASE3_adaptedFromEXIOBASE2.xlsx'),
                'Q_materials'),
            pd.read_excel(
                pkg_resources.resource_stream(__name__,
                                              '/Data/characterisationEXIOBASE3_adaptedFromEXIOBASE2.xlsx'),
                'Q_resources'),
            pd.read_excel(
                pkg_resources.resource_stream(__name__,
                                              '/Data/characterisationEXIOBASE3_adaptedFromEXIOBASE2.xlsx'),
                'Q_factor_inputs')], sort=False).fillna(0)
        self.impact_methods_IO = list(zip(self.C_io.iloc[:, 0], self.C_io.iloc[:, 1], self.C_io.iloc[:, 2],
                                          self.C_io.iloc[:, 3]))
        # reordering + if there is any flows that weren't characterized put a null CF
        self.C_io = self.C_io.reindex(self.flows_of_IO, axis=1).fillna(0)
        self.C_io = scipy.sparse.csr_matrix(self.C_io)

        # reference year allows to determine the inflation rate that will need to be applied
        self.reference_year_IO = int(self.IO_database.meta.description[-4:])

        self.description.append('Exiobase 3 / reference year: '+str(self.reference_year_IO))

        # OPTIONAL FEATURES - MORE ENVIRONMENTAL EXTENSIONS

        if complete_extensions:
            if version_ecoinvent != '3.5':
                raise ValueError("Completed extensions are only available for the hybridization of the 3.5 version "
                                 "of ecoinvent")
            else:
                new_extensions = scipy.sparse.load_npz(pkg_resources.resource_filename(__name__,
                                                                                     '/Data/Completed_extensions_exio3/'
                                                                                     'S_completed.npz'))
                new_flows_names = eval(open(pkg_resources.resource_filename(
                    __name__, '/Data/Completed_extensions_exio3/extension_names.txt'), 'r').read())

                extended_extensions = completing_extensions(pd.DataFrame(self.F_io.todense(),
                                                                            index=self.flows_of_IO,
                                                                            columns=pd.MultiIndex.from_product([
                                                                                  self.regions_of_IO,
                                                                                  self.sectors_of_IO],
                                                                                  names=['region', 'sector'])),
                                                            pd.DataFrame(new_extensions.todense(),
                                                                            index=new_flows_names,
                                                                            columns=pd.MultiIndex.from_product([
                                                                                self.regions_of_IO,
                                                                                self.sectors_of_IO],
                                                                                names=['region', 'sector'])))
                self.extended_flows_names = extended_extensions.index.tolist()
                self.F_io = back_to_sparse(extended_extensions)

                self.description.append('Environmental extensions were completed')

                if not impact_world:
                    self.C_io = scipy.sparse.load_npz(pkg_resources.resource_filename(__name__, '/Data/'
                                                                                     'Completed_extensions_exio3/'
                                                                                     'C_completed_in_CML.npz'))
                    self.extended_impact_names_CML = eval(open(pkg_resources.resource_filename(
                        __name__, '/Data/Completed_extensions_exio3/name_impact_extended_CML.txt'), 'r').read())
                    self.description.append('Classic impact categories were used')
                else:
                    self.C_io = scipy.sparse.load_npz(pkg_resources.resource_filename(
                        __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/Exiobase_not_regionalized.npz'))
                    self.extended_impact_names_IW_exio = eval(open(pkg_resources.resource_filename(
                        __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/not_regionalized_IW+_EXIOBASE.txt'), 'r').read())
                    self.C_f = scipy.sparse.load_npz(pkg_resources.resource_filename(
                        __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/Ecoinvent_not_regionalized.npz'))
                    self.extended_impact_names_IW_eco = eval(open(pkg_resources.resource_filename(
                        __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/not_regionalized_IW+_ecoinvent.txt'), 'r').read())
                    self.description.append('Impact World+ was used')

        else:
            self.description.append('Environmental extensions were not completed')
            if impact_world:
                self.C_io = scipy.sparse.load_npz(pkg_resources.resource_filename(
                    __name__, '/Data/Characterization_matrix_IW+/eco'+version_ecoinvent+'/Exiobase_not_regionalized.npz'))
                self.extended_impact_names_IW_exio = eval(open(pkg_resources.resource_filename(
                    __name__, '/Data/Characterization_matrix_IW+/eco'+version_ecoinvent+'/not_regionalized_IW+_EXIOBASE.txt'), 'r').read())
                self.C_f = scipy.sparse.load_npz(pkg_resources.resource_filename(
                    __name__, '/Data/Characterization_matrix_IW+/eco'+version_ecoinvent+'/Ecoinvent_not_regionalized.npz'))
                self.extended_impact_names_IW_eco = eval(open(pkg_resources.resource_filename(
                    __name__, '/Data/Characterization_matrix_IW+/eco'+version_ecoinvent+'/not_regionalized_IW+_ecoinvent.txt'), 'r').read())
                self.description.append('Impact World+ was used')
                self.F_io = scipy.sparse.load_npz(pkg_resources.resource_filename(
                    __name__, '/Data/Characterization_matrix_IW+/eco'+version_ecoinvent+'/Exiobase_not_regionalized_F.npz'))
                self.extended_flows_names = eval(open(pkg_resources.resource_filename(
                    __name__, '/Data/Characterization_matrix_IW+/eco'+version_ecoinvent+'/normal_flows_names.txt'), 'r').read())
            else:
                self.description.append('Classic impact categories were used')

        # OPTIONAL FEATURES - REGIONALIZATION OF IMPACTS

        if regionalized:
            if version_ecoinvent != '3.5':
                raise ValueError("Regionalization is only available for the hybridization of the 3.5 version "
                                 "of ecoinvent")
            else:
                if not impact_world:
                    print('Only the impact world+ method supports regionalization. '
                          'Please pass the corresponding argument as True if you wish to regionalize')
                else:
                    self.description.append('Regionalized flows/impacts available')
                    if complete_extensions:
                        # load regionalized characterization matrices
                        self.C_io_regio = scipy.sparse.load_npz(pkg_resources.resource_stream(
                            __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/Exiobase_regionalized.npz'))
                        self.regionalized_impact_names_exio = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/regionalized_IW+_exiobase.txt'), 'r').read())
                        self.regionalized_flow_names_exio = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/regionalized_IW+_exiobase_flows.txt'), 'r').read())
                        self.C_f_regio = scipy.sparse.load_npz(pkg_resources.resource_stream(
                            __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/Ecoinvent_regionalized.npz'))
                        self.regionalized_impact_names_eco = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/regionalized_IW+_ecoinvent.txt'), 'r').read())
                        self.regionalized_flow_names_eco = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/completed_extensions/regionalized_IW+_ecoinvent_flows.txt'), 'r').read())
                    else:
                        # load regionalized characterization matrices
                        self.C_io_regio = scipy.sparse.load_npz(pkg_resources.resource_stream(
                            __name__, '/Data/Characterization_matrix_IW+/eco3.5/normal_extensions/Exiobase_regionalized.npz'))
                        self.regionalized_impact_names_exio = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/normal_extensions/regionalized_IW+_exiobase.txt'), 'r').read())
                        self.regionalized_flow_names_exio = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/normal_extensions/regionalized_IW+_exiobase_flows.txt'), 'r').read())
                        self.C_f_regio = scipy.sparse.load_npz(pkg_resources.resource_stream(
                            __name__, '/Data/Characterization_matrix_IW+/eco3.5/normal_extensions/Ecoinvent_regionalized.npz'))
                        self.regionalized_impact_names_eco = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/normal_extensions/regionalized_IW+_ecoinvent.txt'), 'r').read())
                        self.regionalized_flow_names_eco = eval(open(pkg_resources.resource_filename(
                                __name__, '/Data/Characterization_matrix_IW+/eco3.5/normal_extensions/regionalized_IW+_ecoinvent_flows.txt'), 'r').read())
                    # manufacture regionalized environmental extensions matrix (EXIOBASE)
                    F_io = pd.DataFrame(self.F_io.todense(),
                                        self.extended_flows_names,
                                        pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO],
                                                                   names=['region', 'sector']))
                    regionalizable_flows = []
                    for flow in self.regionalized_flow_names_exio:
                        if flow[0] not in regionalizable_flows:
                            regionalizable_flows.append(flow[0])
                    self.F_io_regio = pd.concat([F_io.loc[regionalizable_flows]] * 49)
                    self.F_io_regio.index = pd.MultiIndex.from_product([self.regions_of_IO, regionalizable_flows])
                    self.F_io_regio = self.F_io_regio.swaplevel()
                    F_io_regio_dict = self.F_io_regio.to_dict()

                    for sector in F_io_regio_dict:
                        for flow in F_io_regio_dict[sector].keys():
                            if flow[1] != sector[0]:
                                F_io_regio_dict[sector][flow] = 0

                    self.F_io_regio = pd.DataFrame.from_dict(F_io_regio_dict)
                    self.F_io_regio = self.F_io_regio.reindex(self.regionalized_flow_names_exio)
                    self.F_io_regio = back_to_sparse(self.F_io_regio)
                    del F_io
                    # manufacture environmental flows matrix (ecoinvent)
                    F = pd.DataFrame(self.F_f.todense(), self.STR_f.index, self.PRO_f.index)
                    geography_for_IW = dict(zip(self.PRO_f.index,
                                                [i.split('-')[0].replace('Canada without Quebec', 'CA') for i in
                                                 self.PRO_f.geography]))
                    df_to_pivot = pd.DataFrame()
                    for process in F.columns:
                        regionalizable_flows = [i for i in F.loc[:, process][F.loc[:, process] != 0].index if
                                                i in list(set([i[0] for i in self.regionalized_flow_names_eco]))]
                        for flow in F.loc[regionalizable_flows, process].index:
                            df_to_pivot = pd.concat([df_to_pivot,
                                                     pd.DataFrame([process,
                                                                   geography_for_IW[process],
                                                                   flow,
                                                                   F.loc[regionalizable_flows, process][flow]],
                                                                  index=['process', 'geography', 'flow', 'value']).T],
                                                    sort=False)

                    self.F_f_regio = pd.pivot_table(df_to_pivot, index=['flow', 'geography'], columns='process',
                                                    values='value', aggfunc=sum)
                    self.F_f_regio = self.F_f_regio.fillna(0)
                    self.F_f_regio = pd.concat(
                        [pd.DataFrame(0, self.F_f_regio.index,
                                      [i for i in self.PRO_f.index if i not in self.F_f_regio.columns]),
                         self.F_f_regio], axis=1).reindex(self.PRO_f.index, axis=1)
                    self.F_f_regio = back_to_sparse(self.F_f_regio)
                    del F
                    # remove regionalized characterization from standard characterization matrices
                    self.C_f = pd.DataFrame(self.C_f.todense(), self.extended_impact_names_IW_eco, self.STR_f.MATRIXID)
                    regionalized_flows = set([i[0] for i in self.regionalized_flow_names_eco])
                    for category in self.regionalized_impact_names_eco:
                        # identify flows that have an impact for the given category
                        impacting_flows = self.C_f.loc[category][self.C_f.loc[category] != 0].index.tolist()
                        # if they are characterized with the regionalized matrix, set them to zero
                        self.C_f.loc[
                            category, [i for i in impacting_flows if i in regionalized_flows]] = 0
                    self.C_f = back_to_sparse(self.C_f)
                    self.C_io = pd.DataFrame(self.C_io.todense(), self.extended_impact_names_IW_exio,
                                             self.extended_flows_names)
                    for category in self.regionalized_impact_names_exio:
                        impacting_flows = self.C_io.loc[category][
                            self.C_io.loc[category] != 0].index.tolist()
                        self.C_io.loc[category, [i for i in impacting_flows if i in set(
                            [i[0] for i in self.regionalized_flow_names_exio])]] = 0
                    self.C_io = back_to_sparse(self.C_io)
        else:
            self.description.append('Regionalized flows/impacts unavailable')

        # OPTIONAL FEATURES - CAPITAL GOODS

        if path_to_capitals == '':
            print('No path for the capital folder was provided. Capitals will not be endogenized')
            Kbar = None
            inv_diag_x = None
        else:
            if version_ecoinvent != '3.5':
                raise ValueError("Capital endogenization is only available for the hybridization of the 3.5 version "
                                 "of ecoinvent")
            else:
                Kbar = scipy.io.loadmat(path_to_capitals)['KbarCfc']
                warnings.filterwarnings("ignore", "divide by zero encountered in true_divide")
                inv_diag_x = 1 / self.IO_database.x.values
                inv_diag_x[inv_diag_x == np.inf] = 0
                inv_diag_x = np.diagflat(inv_diag_x)
                self.K_io = Kbar.dot(scipy.sparse.csr_matrix(inv_diag_x))

        del self.LCA_database
        del self.IO_database

        # STAM CATEGORIES

        self.io_categories = ast.literal_eval(pkg_resources.resource_string(__name__,
                                                                            '/Data/eco' + str(
                                                                                version_ecoinvent) + '_exio' + str(
                                                                                version_exiobase) +
                                                                            '/STAM_categories.txt').decode('utf-8'))
        self.categories_same_functionality = ast.literal_eval(
            pkg_resources.resource_string(
                __name__,
                '/Data/eco' + str(version_ecoinvent) + '_exio' + str(
                    version_exiobase) + '/STAM_functional_categories.txt').decode('utf-8'))

        # GEOGRAPHY CONCORDANCE

        self.listcountry = ast.literal_eval(
            pkg_resources.resource_string(__name__,
                '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/countries.txt').decode('utf-8'))

        self.countries_per_regions = ast.literal_eval(
            pkg_resources.resource_string(
                __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
                '/countries_per_regions.txt').decode('utf-8'))

        self.replacements = ast.literal_eval(
            pkg_resources.resource_string(
                __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) +
                '/geography_replacements.txt').decode('utf-8'))
        self.PRO_f['io_geography'] = [self.replacements[i] if i in self.replacements else i for i in self.PRO_f.geography]

        # PRODUCT CONCORDANCE

        concordance_activity = pd.read_excel(pkg_resources.resource_stream(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Product_Concordances.xlsx'),
            'Concordance per activity')
        concordance_product = pd.read_excel(pkg_resources.resource_stream(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Product_Concordances.xlsx'),
            'Concordance per product')
        convert_sector_code = pd.read_excel(pkg_resources.resource_stream(
            __name__,
            '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Product_Concordances.xlsx'),
            'Description_Exiobase')

        concordance_activity = concordance_activity.drop('activityName', axis=1)
        concordance_product = concordance_product.drop('productName', axis=1)
        self.PRO_f = self.PRO_f.merge(concordance_product, 'outer')
        self.PRO_f = self.PRO_f.merge(concordance_activity, how='left', on='activityNameId')
        list_concordance = []
        for index in self.PRO_f.index:
            if type(self.PRO_f.Concordance_y[index]) != float:
                list_concordance.append(self.PRO_f.Concordance_y[index])
            elif type(self.PRO_f.Concordance_y[index]) == float:
                list_concordance.append(self.PRO_f.Concordance_x[index])
        self.PRO_f['Concordance'] = list_concordance
        self.PRO_f = self.PRO_f.drop(['Concordance_x', 'Concordance_y'], axis=1)
        # convert exiobase codes (e.g. p01) to names of the sectors (e.g. 'Paddy rice')
        self.PRO_f = self.PRO_f.merge(convert_sector_code, left_on='Concordance', right_on='EXIO_code', how='left')
        self.PRO_f = self.PRO_f.drop(['Concordance', 'EXIO_code'], axis=1)
        self.PRO_f.index = self.PRO_f.activityId + '_' + self.PRO_f.productId
        # We want the indexes of PRO_f and A_ff in the same order
        self.PRO_f = self.PRO_f.reindex(self.processes_in_order)

        # LOADING THE FILTER

        self.list_to_hyb = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Hybridized', index_col=0).index.tolist()
        self.listmarket = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Market', index_col=0).index.tolist()
        self.listnottransacted = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Not commercialized', index_col=0).index.tolist()
        self.listguillotine = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Poor quality', index_col=0).index.tolist()
        self.dummyprocesses = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Empty processes', index_col=0).index.tolist()
        self.null_price = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'No price', index_col=0).index.tolist()
        self.list_uncovered_geographies = pd.read_excel(pkg_resources.resource_stream(
            __name__, '/Data/eco' + str(version_ecoinvent) + '_exio' + str(version_exiobase) + '/Filter.xlsx'),
            'Uncovered geography', index_col=0).index.tolist()
        self.list_not_to_hyb = (
                self.listmarket + self.listnottransacted + self.listguillotine + self.dummyprocesses
                + self.null_price + self.list_uncovered_geographies)

        return LCAIO(PRO_f=self.PRO_f, A_ff=self.A_ff, A_io=self.A_io, F_f=self.F_f, F_io=self.F_io, X_io=self.X_io,
                     y_io=self.y_io, C_f=self.C_f, C_io=self.C_io, STR_f=self.STR_f, IMP=self.IMP,
                     listcountry=self.listcountry, K_io=self.K_io,
                     countries_per_regions=self.countries_per_regions, reference_year_IO=self.reference_year_IO,
                     number_of_countries_IO=self.number_of_countries_IO, number_of_RoW_IO=self.number_of_RoW_IO,
                     number_of_products_IO=self.number_of_products_IO, list_to_hyb=self.list_to_hyb,
                     list_not_to_hyb=self.list_not_to_hyb, listmarket=self.listmarket,
                     dummyprocesses=self.dummyprocesses, listnottransacted=self.listnottransacted,
                     null_price=self.null_price, listguillotine=self.listguillotine,
                     list_uncovered_geographies=self.list_uncovered_geographies, io_categories=self.io_categories,
                     categories_same_functionality=self.categories_same_functionality,
                     lca_database_name_and_version=self.lca_database_name_and_version,
                     io_database_name_and_version=self.io_database_name_and_version,
                     processes_in_order=self.processes_in_order, sectors_of_IO=self.sectors_of_IO,
                     regions_of_IO=self.regions_of_IO, Y_categories=self.Y_categories, flows_of_IO=self.flows_of_IO,
                     impact_methods_IO=self.impact_methods_IO, extended_flows_names=self.extended_flows_names,
                     extended_impact_names_CML=self.extended_impact_names_CML,
                     extended_impact_names_IW_exio=self.extended_impact_names_IW_exio,
                     extended_impact_names_IW_eco=self.extended_impact_names_IW_eco, C_io_regio=self.C_io_regio,
                     C_f_regio=self.C_f_regio, F_io_regio=self.F_io_regio, F_f_regio=self.F_f_regio,
                     regionalized_impact_names_exio=self.regionalized_impact_names_exio,
                     regionalized_impact_names_eco=self.regionalized_impact_names_eco,
                     regionalized_flow_names_exio=self.regionalized_flow_names_exio,
                     regionalized_flow_names_eco=self.regionalized_flow_names_eco,  description=self.description)


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
        * correct_inconsistencies()
        * low_production_volume_processes()
        * extend_inventory()
        * extract_scaling_vector_technosphere()
        * extract_scaling_vector_biosphere
        * extract_flow_amounts_technosphere()
        * extract_flow_amounts_biosphere()
        * apply_scaling_without_prices()
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
        self.IMP = pd.DataFrame()
        self.X_io = pd.DataFrame()
        self.C_io_regio = pd.DataFrame()
        self.C_f_regio = pd.DataFrame()
        self.F_io_regio = pd.DataFrame()
        self.F_f_regio = pd.DataFrame()

        self.K_io = pd.DataFrame()
        self.K_io_f = pd.DataFrame()

        self.io_categories = defaultdict(list)
        self.categories_same_functionality = []
        self.processes_in_order = []
        self.sectors_of_IO = []
        self.regions_of_IO = []
        self.Y_categories = []
        self.flows_of_IO = []
        self.impact_methods_IO = []
        self.extended_flows_names = []
        self.extended_impact_names_CML = []
        self.extended_impact_names_IW_exio = []
        self.extended_impact_names_IW_eco = []
        self.regionalized_impact_names_exio = []
        self.regionalized_impact_names_eco = []
        self.regionalized_flow_names_exio = []
        self.regionalized_flow_names_eco = []

        self.listcountry = []
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
        self.hybridized_processes = []

        self.A_ff_processed = pd.DataFrame()
        self.total_prod_country = pd.DataFrame()
        self.total_prod_region = pd.DataFrame()
        self.total_prod_RoW = pd.DataFrame()
        self.dictRoW = {}
        self.STAM_table = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/STAM_table.xlsx'), index_col=0 )
        self.patching_exiobase = pd.read_excel(pkg_resources.resource_stream(__name__, '/Data/'
                                                                                       'Exiobase_patchwork.xlsx'), index_col=0)
        self.H = pd.DataFrame()
        self.G = pd.DataFrame()
        self.A_io_f_uncorrected = pd.DataFrame()
        self.aggregated_A_io = pd.DataFrame()
        self.aggregated_F_io = pd.DataFrame()
        self.add_on_H_scaled_vector = pd.DataFrame()

        self.double_counting = ''
        self.description = []
        self.basic_A_io_f_uncorrected = pd.DataFrame()

        allowed_keys = list(self.__dict__.keys())
        self.__dict__.update((key, value) for key, value in kwargs.items() if key in allowed_keys)

    # ----------------------------CORE METHOD-------------------------------------

    def hybridize(self, method_double_counting, capitals=False):
        """ Hybridize an LCA database with an IO database

        self.A_io_f_uncorrected is calculated following the equation (1) of the paper [insert doi]
        self.A_io_f is calculated following the equation (2) of the paper [insert doi]

        Args:
        -----
            method_double_counting  : method to correct double counting with (='binary' or ='STAM')

        Returns:
        -------
            The updated self.A_io_f matrix

        """

        if not method_double_counting:
            raise Exception('Please enter a method to correct double counting (i.e. binary or STAM)')

        if not capitals:
            self.capitals = False
            self.description.append('Capitals were not endogenized')
        if capitals:
            self.description.append('Capitals were endogenized')
            self.capitals = True

        self.identify_rows()
        self.update_prices_electricity()
        self.calc_productions()
        self.low_production_volume_processes()
        self.extend_inventory()

        # ---- CONVERSION PART ------

        self.H = pd.DataFrame(0, index=self.sectors_of_IO, columns=self.PRO_f.index, dtype='int64')
        for sector in self.H.index:
            self.H.loc[sector, self.H.columns.intersection(self.PRO_f[self.PRO_f.ProductTypeName == sector].index)] = 1
        self.H = self.H.append([self.H] * (self.number_of_countries_IO + self.number_of_RoW_IO - 1))
        self.H.index = pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO], names=['region', 'sector'])

        # translate the geography concordance txt files into matrices
        matrix_countries_per_region = pd.DataFrame(0, index=self.listcountry,
                                                   columns=self.listcountry +
                                                           list(self.countries_per_regions) +
                                                           list(self.dictRoW.keys()))
        for country in self.listcountry:
            matrix_countries_per_region.loc[country, country] = 1
        for region in self.countries_per_regions.keys():
            matrix_countries_per_region.loc[self.countries_per_regions[region], region] = 1
        for RoW in self.dictRoW.keys():
            matrix_countries_per_region.loc[self.dictRoW[RoW], RoW] = 1

        # identify which region corresponds to which process
        region_covered_per_process = pd.DataFrame(0, index=self.listcountry +
                                                          list(self.countries_per_regions) +
                                                          list(self.dictRoW.keys()),
                                                  columns=self.PRO_f.index)
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
                                  self.sectors_of_IO[commodity_position]))
        concordance_geography = pd.concat([matrix_countries_per_region] * self.number_of_products_IO).sort_index()
        concordance_geography.index = listindex
        self.H.index = self.H.index.values
        concordance_geography = concordance_geography.reindex(self.H.index)
        concordance_geography.index = pd.MultiIndex.from_product(
            [self.regions_of_IO, self.sectors_of_IO], names=['region', 'sector'])
        self.H.index = pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO], names=['region', 'sector'])

        # introduce production volumes in the mix
        concordance_geography_with_production_volumes = concordance_geography.multiply(
            pd.concat([self.total_prod_country] *
                      len(region_covered_per_process),
                      axis=1).values)
        # pivot the total_prod_country multiindex serie into a single index dataframe
        pivoting = pd.pivot_table(data=self.total_prod_country, values=self.total_prod_country,
                                  columns=self.total_prod_country.index.get_level_values('region'),
                                  index=self.total_prod_country.index.get_level_values('sector'))
        # some more indexing
        pivoting = pivoting.reindex(
            self.total_prod_country.index.get_level_values('sector')[:self.number_of_products_IO])
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

        # product concordance matrix filtered
        H_for_hyb = self.H.copy()
        H_for_hyb.loc[:, self.list_not_to_hyb] = 0

        self.A_io_f_uncorrected = pd.DataFrame(self.A_io.todense(),
                                               index=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector']),
                                               columns=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector'])).dot(
            H_for_hyb * inflation * Geo) * self.PRO_f.price

        if capitals:
            self.K_io_f_uncorrected = pd.DataFrame(self.K_io.todense(), index=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector']), columns=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector'])).dot(
                H_for_hyb * inflation * Geo) * self.PRO_f.price

            # Since we endogenized capitals we need to remove them from both final demand and factors of production
            self.y_io = pd.DataFrame(self.y_io.todense(),
                                     pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO]),
                                     pd.MultiIndex.from_product([self.regions_of_IO, self.Y_categories]))
            self.y_io.loc[:, [i for i in self.y_io.columns if i[1] == 'Gross fixed capital formation']] = 0
            self.y_io = back_to_sparse(self.y_io)
            if self.extended_flows_names:
                self.F_io = pd.DataFrame(self.F_io.todense(), self.extended_flows_names,
                                         pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO]))
            else:
                self.F_io = pd.DataFrame(self.F_io.todense(), self.flows_of_IO,
                                         pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO]))
            self.F_io.loc['Operating surplus: Consumption of fixed capital', :] = 0
            self.F_io = back_to_sparse(self.F_io)

        # ---- HYBRIDIZATION WITHOUT PRICES ------
        self.apply_scaling_without_prices(self.capitals)
        del self.aggregated_A_io
        del self.aggregated_F_io
        self.A_io_f_uncorrected += pd.DataFrame(self.A_io.todense(), index=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector']), columns=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector'])).dot(
            self.add_on_H_scaled_vector * inflation * Geo)
        self.A_io_f_uncorrected = back_to_sparse(self.A_io_f_uncorrected)

        if capitals:
            self.K_io_f_uncorrected = self.K_io_f_uncorrected + pd.DataFrame(self.K_io.todense(),
                                                                             index=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector']), columns=pd.MultiIndex.from_product(
                                                   [self.regions_of_IO, self.sectors_of_IO],
                                                   names=['region', 'sector'])).dot(
                self.add_on_H_scaled_vector * inflation * Geo)
            self.K_io_f_uncorrected = back_to_sparse(self.K_io_f_uncorrected)

        # ------ DOUBLE COUNTING PART -------

        if method_double_counting == 'binary':
            lambda_filter_matrix = self.H.dot(pd.DataFrame(self.A_ff_processed.todense(),
                                              self.PRO_f.index, self.PRO_f.index))
            lambda_filter_matrix = lambda_filter_matrix.mask(lambda_filter_matrix > 0)
            lambda_filter_matrix[lambda_filter_matrix == 0] = 1
            lambda_filter_matrix = lambda_filter_matrix.fillna(0)
            # if capitals_method == 'LCA' or capitals_method == 'IO':
            #     # double counting of capitals must be corrected in K_io_f and not A_io_f
            #     lambda_filter_matrix.loc[[i for i in lambda_filter_matrix.index
            #                               if i[1] in self.list_of_capital_sectors]] = 1
            lambda_filter_matrix = back_to_sparse(lambda_filter_matrix)

            self.A_io_f = self.A_io_f_uncorrected.multiply(lambda_filter_matrix)
            if capitals:
                self.K_io_f = self.K_io_f_uncorrected.multiply(lambda_filter_matrix)

            self.double_counting = 'binary'
            self.description.append('The binary method was used to correct for double counting')

        elif method_double_counting == 'STAM':
            lambda_filter_matrix = self.H.dot(pd.DataFrame(self.A_ff_processed.todense(),
                                              self.PRO_f.index, self.PRO_f.index))
            lambda_filter_matrix = lambda_filter_matrix.mask(lambda_filter_matrix > 0)
            lambda_filter_matrix[lambda_filter_matrix == 0] = 1
            lambda_filter_matrix = lambda_filter_matrix.fillna(0)
            lambda_filter_matrix = back_to_sparse(lambda_filter_matrix)

            self.G = pd.DataFrame(0, index=self.sectors_of_IO, columns=self.STAM_table.columns)
            for columns in self.G.columns:
                self.G.loc[[i for i in self.G.index if i in self.io_categories[columns]], columns] = 1
            self.G = self.G.append([self.G] * (self.number_of_countries_IO + self.number_of_RoW_IO - 1))
            self.G.index = pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO],
                                                      names=['region', 'sector'])

            gamma_filter_matrix = self.G.dot((self.STAM_table.mul(self.patching_exiobase)).dot(self.G.transpose().
                                                                                               dot(self.H)))
            gamma_filter_matrix[gamma_filter_matrix == self.number_of_countries_IO + self.number_of_RoW_IO] = 1
            gamma_filter_matrix = back_to_sparse(gamma_filter_matrix)

            phi_filter_matrix = pd.DataFrame(1, index=self.G.index, columns=self.PRO_f.index)
            categories_used_by_processes = self.G.transpose().dot(self.H.dot(
                pd.DataFrame(self.A_ff_processed.todense(), self.PRO_f.index, self.PRO_f.index)))
            del self.G

            for category in self.categories_same_functionality:
                list_category_to_zero = [i for i in categories_used_by_processes.columns if
                                         categories_used_by_processes.loc[category, i] != 0]
                phi_filter_matrix.loc[[i for i in phi_filter_matrix.index if
                                       i[1] in self.io_categories[category]], list_category_to_zero] = 0
            del categories_used_by_processes
            phi_filter_matrix = back_to_sparse(phi_filter_matrix)

            self.A_io_f = phi_filter_matrix.multiply(
                gamma_filter_matrix.multiply(lambda_filter_matrix.multiply(self.A_io_f_uncorrected)))

            if capitals:
                self.K_io_f = phi_filter_matrix.multiply(
                    gamma_filter_matrix.multiply(lambda_filter_matrix.multiply(self.K_io_f_uncorrected)))

            self.double_counting = 'STAM'
            self.description.append('STAM was used to correct for double counting')

    # ---------------------PREPARATIONS FOR THE HYBRIDIZATION----------------------

    def identify_rows(self):
        """ Method to identify the various unique Rest of the World (RoW) regions of the LCA database

        Returns:
        --------
            The updated self.dictRoW in which unique RoW regions are identified in terms of countries they include

        """

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
        self.PRO_f = self.PRO_f.reindex(self.processes_in_order)
        self.PRO_f.io_geography.update(self.PRO_f.RoW_geography[self.PRO_f.io_geography == 'RoW'])
        self.PRO_f = self.PRO_f.drop('RoW_geography', axis=1)
        # might be some RoW or empty lists left in PRO_f
        self.PRO_f.io_geography[self.PRO_f.io_geography == 'RoW'] = 'GLO'
        self.PRO_f.io_geography.loc[[i for i in self.PRO_f.index if type(self.PRO_f.io_geography[i]) == list]] = 'GLO'

    def update_prices_electricity(self):
        """ Method specially for ecoinvent and exiobase, in which electricity prices from ecoinvent are replaced by
            price extracted from exiobase

         Returns:
         -------
            The updated price column of the self.PRO_f matrix

         """

        electricity_price = pd.read_excel(pkg_resources.resource_stream(__name__,
                                                                        '/Data/Regionalized_electricity_prices.xlsx'),
                                                                        index_col=0)

        electricity_processes = self.PRO_f.price.loc[
            [i for i in self.PRO_f.index if
             (0.0977 == self.PRO_f.price[i] and 'electricity' in self.PRO_f.productName[i]) or
             (0.107 == self.PRO_f.price[i] and 'electricity' in self.PRO_f.productName[i])]]

        if electricity_processes.empty:
            print('Empty!')
            return
        merged = self.PRO_f.loc[electricity_processes.index.values, ['price', 'io_geography']].merge(
            electricity_price, left_on=['io_geography'], right_on=electricity_price.index, how='left')
        merged.prices[merged.prices.isnull()] = merged.price
        merged.index = electricity_processes.index.values
        merged = merged.drop(['price', 'io_geography'], axis=1)
        self.PRO_f.price.update(merged.prices)

    def calc_productions(self):
        """ Calculates the different total productions for either countries, regions or RoWs

        Returns:
        -------
            The updated self.total_prod_country, self.total_prod_region and self.total_prod_RoW dataframe

        """

        # the user needs to determine the total demand before being able to calculate productions
        listdrop = []

        absent_countries = {}
        for i in range(0, len(list(self.countries_per_regions.values()))):
            absent_country = [item for item in self.listcountry if
                              item not in list(self.countries_per_regions.values())[i]]
            absent_countries[list(self.countries_per_regions.keys())[i]] = absent_country

        self.total_prod_country = self.total_prod_country = pd.DataFrame(
            self.X_io.todense(), index=pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO],
                                                                  names=['region', 'sector']), columns=['production'])

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

    def correct_inconsistencies(self):
        """ While introducing the completed extensions of EXIOBASE based on USEEIO, some missing emissions are added
        to sectors which in reality produce no products. For these sectors, the added emissions are put back to zero."""

        null_production_sectors = [i for i in pd.MultiIndex.from_product(
            [self.regions_of_IO, self.sectors_of_IO]) if self.total_prod_country.loc[i][0] == 0]

        if 'Environmental extensions were completed' in self.description:
            self.F_io = pd.DataFrame(self.F_io.todense(), self.extended_flows_names,
                                     pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO]))
            self.F_io = self.F_io.drop(self.F_io.loc[:, null_production_sectors].columns, axis=1)
            self.F_io = pd.concat(
                [self.F_io, pd.DataFrame(0, index=self.F_io.index, columns=null_production_sectors)], axis=1)
            self.F_io = self.F_io.reindex(pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO]), axis=1)
            self.F_io = back_to_sparse(self.F_io)

        if 'Regionalized flows/impacts available' in self.description:
            self.F_io_regio = pd.DataFrame(self.F_io_regio.todense(), self.regionalized_flow_names_exio,
                                           pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO],
                                                                      names=['region', 'sector']))
            self.F_io_regio.loc[:, null_production_sectors] = 0
            self.F_io_regio = back_to_sparse(self.F_io_regio)

    def low_production_volume_processes(self):

        list_low_prod_sectors = self.total_prod_country[self.total_prod_country < 10].dropna().index.tolist()

        df = self.PRO_f.loc[[i for i in self.PRO_f.index if (self.PRO_f.io_geography[i], self.PRO_f.ProductTypeName[i])
                             in list_low_prod_sectors and i in self.list_to_hyb]]
        dict_ = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Classing_countries.txt').decode('utf-8'))
        for process in self.PRO_f.index:
            if process in df.index:
                if dict_[self.PRO_f.io_geography[process]] == 'RER':
                    self.PRO_f.io_geography[process] = dict_[self.PRO_f.io_geography[process]]

    def extend_inventory(self):
        """ Method creating a new technology matrix for the LCA database in which the inputs of processes not to
        hybridize are added to the unit process inventories of each LCA process to hybridize

        Returns:
        -------
            The updated self.A_ff_processed matrix

        """

        # matrix of non-hybridized processes
        ANT = self.A_ff.copy()
        ANT = pd.DataFrame(ANT.todense(), self.PRO_f.index, self.PRO_f.index)
        ANT.loc[self.list_to_hyb] = 0
        ANT.fillna(0,inplace=True)

        # matrix of production processes
        Amarket = self.A_ff.copy()
        Amarket = pd.DataFrame(Amarket.todense(), self.PRO_f.index, self.PRO_f.index)
        Amarket.loc[self.listmarket] = 0
        Amarket.fillna(0,inplace=True)

        self.A_ff_processed = Amarket.values.dot(np.linalg.inv(np.eye(len(ANT.values)) - ANT.values))

        # the inversion brings a lot of noise which falsely contributes to double counting which is why it is removed
        self.A_ff_processed[self.A_ff_processed < 10 ** -16] = 0
        self.A_ff_processed = back_to_sparse(self.A_ff_processed)

    def extract_scaling_vector_technosphere(self, list_add_on_to_hyb, scaling_flow, sector_of_scaling_flow,
                                            sector_of_add_ons):
        """This functions extracts the scaling vector used to bypass the poor price data quality of some processes
        of ecoinvent
        Args:
        ----
            list_add_on_to_hyb: a list of uuids of processes for which scaling vectors are needed
            scaling_flow: the name of the flow on which to scale ('concrete','steel','kerosene','aluminium',
                         'electricity')
            sector_of_scaling_flow: the economic sector of the scaling flow
            sector_of_add_ons: the economic sector of the processes that are part of the add_on

        """

        add_on_H_for_hyb = self.H.copy()
        add_on_H_for_hyb.loc[:, [i for i in add_on_H_for_hyb.columns if i not in list_add_on_to_hyb]] = 0

        dff = self.extract_flow_amounts_technosphere(list_add_on_to_hyb, scaling_flow)

        scaling_vector = self.PRO_f.price.copy()
        scaling_vector.loc[:] = 0
        if sector_of_scaling_flow == 'Energy':
            list_energies = ['Electricity by Geothermal', 'Electricity by biomass and waste', 'Electricity by coal',
                             'Electricity by gas', 'Electricity by hydro', 'Electricity by nuclear',
                             'Electricity by petroleum and other oil derivatives', 'Electricity by solar photovoltaic',
                             'Electricity by solar thermal', 'Electricity by tide, wave, ocean', 'Electricity by wind',
                             'Electricity nec', 'Gas/Diesel Oil']
            scaling_vector.loc[dff.index] = (dff / (self.aggregated_A_io.loc[list_energies, sector_of_add_ons].sum() /
                                                    (self.number_of_countries_IO+self.number_of_RoW_IO))).iloc[:, 0]
        elif sector_of_scaling_flow == 'Electricity':
            list_energies = ['Electricity by Geothermal', 'Electricity by biomass and waste', 'Electricity by coal',
                             'Electricity by gas', 'Electricity by hydro', 'Electricity by nuclear',
                             'Electricity by petroleum and other oil derivatives', 'Electricity by solar photovoltaic',
                             'Electricity by solar thermal', 'Electricity by tide, wave, ocean', 'Electricity by wind',
                             'Electricity nec']
            scaling_vector.loc[dff.index] = (dff / (self.aggregated_A_io.loc[list_energies, sector_of_add_ons].sum() /
                                                    (self.number_of_countries_IO+self.number_of_RoW_IO))).iloc[:, 0]
        else:
            scaling_vector.loc[dff.index] = (dff / (self.aggregated_A_io.loc[sector_of_scaling_flow,
                                                                             sector_of_add_ons] /
                                                    (self.number_of_countries_IO+self.number_of_RoW_IO))).iloc[:, 0]

        self.add_on_H_scaled_vector = pd.concat([self.add_on_H_scaled_vector.drop(list_add_on_to_hyb, axis=1),
                                                 (add_on_H_for_hyb.multiply(scaling_vector)).loc[:,
                                                 list_add_on_to_hyb]], axis=1).reindex(columns=self.H.columns)

    def extract_scaling_vector_biosphere(self, what_is_treated, treatment):
        """This functions extracts the scaling vector used to bypass the poor price data quality of some processes
        of ecoinvent
        Args:
        ----
            list_add_on_to_hyb: a list of uuids of processes for which scaling vectors are needed
            scaling_flow: the name of the flow on which to scale ('concrete','steel','kerosene','aluminium',
                         'electricity')
            sector_of_scaling_flow: the economic sector of the scaling flow
            sector_of_add_ons: the economic sector of the processes that are part of the add_on

        """
        add_on_H_for_hyb = self.H.copy()
        list_add_on_to_hyb = [i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                              == str(what_is_treated) + ' waste for treatment: ' + str(treatment) and (
                                      i in self.listguillotine or self.null_price) and
                              'market for' not in self.PRO_f.activityName[i] and
                              i not in self.list_to_hyb and
                              'open' not in self.PRO_f.activityName[i]]
        add_on_H_for_hyb.loc[:, [i for i in add_on_H_for_hyb.columns if i not in list_add_on_to_hyb]] = 0

        biogenic_sectors = ['Other waste for treatment: waste water treatment',
                            'Food waste for treatment: waste water treatment',
                            'Paper waste for treatment: landfill',
                            'Food waste for treatment: landfill',
                            'Wood waste for treatment: landfill',
                            'Food waste for treatment: composting and land application',
                            # 'Wood waste for treatment: incineration',
                            # 'Food waste for treatment: incineration',
                            # 'Paper waste for treatment: incineration'
                            ]

        if str(what_is_treated) + ' waste for treatment: ' + str(treatment) in biogenic_sectors:
            dff = self.extract_flow_amounts_biosphere(list_add_on_to_hyb, True)
            if self.extended_flows_names:
                list_CO2_extensions = [i for i in self.extended_flows_names if 'CO2' in i and 'biogenic' in i]
            else:
                list_CO2_extensions = [i for i in self.flows_of_IO if 'CO2' in i and 'biogenic' in i]
        else:
            dff = self.extract_flow_amounts_biosphere(list_add_on_to_hyb, False)
            if self.extended_flows_names:
                list_CO2_extensions = [i for i in self.extended_flows_names if 'CO2' in i and 'biogenic' not in i]
            else:
                list_CO2_extensions = [i for i in self.flows_of_IO if 'CO2' in i and 'biogenic' not in i]

        scaling_vector = self.PRO_f.price.copy()
        scaling_vector.loc[:] = 0
        if what_is_treated == 'Sewage sludge':
            scaling_vector.loc[dff.index] = (dff / (self.aggregated_F_io.loc[list_CO2_extensions, str(what_is_treated) +
                                                                             ' for treatment: ' + str(treatment)].
                                                    sum() /(self.number_of_countries_IO+self.number_of_RoW_IO))).iloc[:, 0]
        else:
            scaling_vector.loc[dff.index] = (dff / (self.aggregated_F_io.loc[list_CO2_extensions, str(what_is_treated) +
                                                                             ' waste for treatment: ' + str(treatment)].
                                                    sum() /(self.number_of_countries_IO+self.number_of_RoW_IO))).iloc[:, 0]

        self.add_on_H_scaled_vector = pd.concat([self.add_on_H_scaled_vector.drop(list_add_on_to_hyb, axis=1),
                                                 (add_on_H_for_hyb.multiply(scaling_vector)).loc[:,
                                                 list_add_on_to_hyb]], axis=1).reindex(columns=self.H.columns)

    def extract_flow_amounts_technosphere(self, list_of_uuids, flow):
        """This function extracts the amount (in euros) necessary for the determination of a scaling vector
        from the list of uuids passed as an argument. Only flows from the technosphere are accepted in this function.

        Args:
        ----
            list_of_uuids: a list of UUIDs from which to extract the scaling vector values
            flow: the name of the flow on which to scale ('concrete','steel','kerosene','aluminium','energy')

        Returns:
        -------
            a dataframe of the values of the desired flow (concrete, steel, etc.) per UUID
        """
        if flow == 'concrete':
            keyword = [i for i in self.PRO_f.index if
                       ('concrete' in self.PRO_f.productName[i] or 'cement, unspecified' in self.PRO_f.productName[i])]
            second_keyword = [i for i in self.PRO_f.index if ('building, hall') in self.PRO_f.productName[i]]
            return pd.DataFrame(
                self.A_ff.loc[keyword, list_of_uuids].mul(self.PRO_f.price.loc[keyword], axis=0).sum(axis=0) + (
                        self.A_ff.loc[second_keyword, list_of_uuids] * 0.3 * 185).sum(axis=0))
        elif flow == 'steel':
            keyword = [i for i in self.PRO_f.index if
                       ('steel' in self.PRO_f.productName[i] or 'cast iron' in self.PRO_f.productName[i])]
            return pd.DataFrame(
                self.A_ff.loc[keyword, list_of_uuids].mul(self.PRO_f.price.loc[keyword], axis=0).sum(axis=0))
        elif flow == 'energy':
            keyword = [i for i in self.PRO_f.index if
                       ('electricity' in self.PRO_f.productName[i] or 'diesel' in self.PRO_f.productName[i])]
            return pd.DataFrame(
                self.A_ff.loc[keyword, list_of_uuids].mul(self.PRO_f.price.loc[keyword], axis=0).sum(axis=0))
        elif flow == 'electricity':
            keyword = [i for i in self.PRO_f.index if 'electricity' in self.PRO_f.productName[i]]
            return pd.DataFrame(
                self.A_ff.loc[keyword, list_of_uuids].mul(self.PRO_f.price.loc[keyword], axis=0).sum(axis=0))
        else:
            keyword = [i for i in self.PRO_f.index if flow in self.PRO_f.productName[i]]
            return pd.DataFrame(
                self.A_ff.loc[keyword, list_of_uuids].mul(self.PRO_f.price.loc[keyword], axis=0).sum(axis=0))

    def extract_flow_amounts_biosphere(self, list_of_uuids, biogenic=True):
        """This function extracts the amount (in euros) necessary for the determination of a scaling vector
        from the list of uuids passed as an argument. Only flows from the biosphere are accepted in this function.
        By default, the flow will be CO2 (for now).

        Args:
        ----
            list_of_uuids   : a list of UUIDs from which to extract the scaling vector values
            biogenic        : boolean, if biogenic carbon is to be added to the total amount of carbon or not

        Returns:
        -------
            a dataframe of the values of the desired flow (concrete, steel, etc.) per UUID
        """
        if biogenic:
            keyword = [i for i in self.F_f.index if
                       ('Carbon dioxide' in self.STR_f.FULLNAME[i] and 'non-fossil' in self.STR_f.FULLNAME[i])]
        else:
            keyword = [i for i in self.F_f.index if
                       ('Carbon dioxide' in self.STR_f.FULLNAME[i] and 'non-fossil' not in self.STR_f.FULLNAME[i])]
        return pd.DataFrame(self.F_f.loc[keyword, list_of_uuids].sum())

    def apply_scaling_without_prices(self, capitals):
        # The following lines of code are introducing add-ons to enable hybridization bypassing bad price quality
        self.add_on_H_scaled_vector = pd.DataFrame(0, self.H.index, self.H.columns)

        if capitals:
            # keep A+K cause we are looking at all the concrete (so renovation + construction)
            self.aggregated_A_io = (
                    pd.DataFrame(self.A_io.todense(), self.H.index, self.H.index) +
                    pd.DataFrame(self.K_io.todense(), self.H.index, self.H.index)).groupby('sector', axis=1).sum()
        else:
            self.aggregated_A_io = pd.DataFrame(
                self.A_io.todense(), self.H.index, self.H.index).groupby('sector', axis=1).sum()
        self.aggregated_A_io = self.aggregated_A_io.groupby('sector', axis=0).sum()
        # identify if environmental flows got extended or not in database_loader()
        if self.extended_flows_names:
            self.aggregated_F_io = pd.DataFrame(self.F_io.todense(), self.extended_flows_names, self.H.index).groupby(
                level=1, axis=1).sum()
        else:
            self.aggregated_F_io = pd.DataFrame(self.F_io.todense(), self.flows_of_IO, self.H.index).groupby(
                level=1, axis=1).sum()

        self.A_ff = pd.DataFrame(self.A_ff.todense(), self.PRO_f.index, self.PRO_f.index)
        self.F_f = pd.DataFrame(self.F_f.todense(), self.STR_f.index, self.PRO_f.index)

        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Construction work' and i in self.listguillotine],
                                                 'concrete',
                                                 'Cement, lime and plaster',
                                                 'Construction work')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Construction work' and i in self.listguillotine])

        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Other transport equipment' and i in self.listguillotine
                                                  and 'aircraft' not in self.PRO_f.productName[i]],
                                                 'steel',
                                                 'Basic iron and steel and of ferro-alloys and first products thereof',
                                                 'Other transport equipment')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Other transport equipment' and i in self.listguillotine
                                                  and 'aircraft' not in self.PRO_f.productName[i]])

        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Other transport equipment' and i in self.listguillotine and
                                                  'aircraft' in self.PRO_f.productName[i]],
                                                 'aluminium',
                                                 'Aluminium and aluminium products',
                                                 'Other transport equipment')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Other transport equipment' and i in self.listguillotine and
                                                  'aircraft' in self.PRO_f.productName[i]])

        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Machinery and equipment n.e.c.' and i in self.listguillotine],
                                                 'steel',
                                                 'Basic iron and steel and of ferro-alloys and first products thereof',
                                                 'Machinery and equipment n.e.c.')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                            'Machinery and equipment n.e.c.' and i in self.listguillotine])

        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Air transport services' and i in self.listguillotine],
                                                 'kerosene',
                                                 'Kerosene Type Jet Fuel',
                                                 'Air transport services')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Air transport services' and i in self.listguillotine])

        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Railway transportation services' and i in self.listguillotine],
                                                 'energy',
                                                 'Energy',
                                                 'Railway transportation services')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i] ==
                                                  'Railway transportation services' and i in self.listguillotine])

        self.extract_scaling_vector_biosphere('Food', 'incineration')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Food waste for treatment: incineration' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Paper', 'incineration')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Paper waste for treatment: incineration' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Plastic', 'incineration')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Plastic waste for treatment: incineration' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Intert/metal', 'incineration')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Intert/metal waste for treatment: incineration' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Textiles', 'incineration')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Textiles waste for treatment: incineration' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Wood', 'incineration')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Wood waste for treatment: incineration' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Oil/hazardous', 'incineration')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Oil/hazardous waste for treatment: incineration' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])

        self.extract_scaling_vector_biosphere('Plastic', 'landfill')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Plastic waste for treatment: landfill' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Inert/metal/hazardous', 'landfill')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Inert/metal/hazardous waste for treatment: landfill' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Wood', 'landfill')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Wood waste for treatment: landfill' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])

        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                                  == 'Food waste for treatment: waste water treatment' and (
                                                          i in self.listguillotine or self.null_price) and
                                                  'market for' not in self.PRO_f.activityName[i] and
                                                  i not in self.list_to_hyb and
                                                  'open' not in self.PRO_f.activityName[i]],
                                                 'electricity',
                                                 'Electricity',
                                                 'Food waste for treatment: waste water treatment')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Food waste for treatment: waste water treatment' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_technosphere([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                                  == 'Other waste for treatment: waste water treatment' and (
                                                          i in self.listguillotine or self.null_price) and
                                                  'market for' not in self.PRO_f.activityName[i] and
                                                  i not in self.list_to_hyb and
                                                  'open' not in self.PRO_f.activityName[i]],
                                                 'electricity',
                                                 'Electricity',
                                                 'Other waste for treatment: waste water treatment')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Other waste for treatment: waste water treatment' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])

        self.extract_scaling_vector_biosphere('Food', 'biogasification and land application')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Food waste for treatment: biogasification and land application' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Sewage sludge', 'biogasification and land application')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Sewage sludge waste for treatment: biogasification and land application' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])

        self.extract_scaling_vector_biosphere('Food', 'composting and land application')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Food waste for treatment: composting and land application' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])
        self.extract_scaling_vector_biosphere('Paper and wood', 'composting and land application')
        self.hybridized_processes.append([i for i in self.H.columns if self.PRO_f.ProductTypeName[i]
                                          == 'Paper and wood waste for treatment: composting and land application' and (
                                                  i in self.listguillotine or self.null_price) and
                                          'market for' not in self.PRO_f.activityName[i] and
                                          i not in self.list_to_hyb and 'open' not in self.PRO_f.activityName[i]])

        self.hybridized_processes = [item for sublist in self.hybridized_processes for item in sublist]
        self.hybridized_processes += self.list_to_hyb

        self.A_ff = back_to_sparse(self.A_ff)
        self.F_f = back_to_sparse(self.F_f)

    # -------------------------- EXPORT RESULTS -----------------------------------

    def save_system(self, format='pickle'):
        """ Export the hybridized database to dataframe via pickle or csv

        Args:
        -----
                format: 'pickle' by default, find it in the Databases folder of the module
        Returns:
        --------
            A pickle of the hybrid system exported to the corresponding folder in the pylcaio package

        """

        if 'Capitals were endogenized' in self.description:
            if ('Impact World+ was used' in self.description and
                    'Regionalized flows/impacts available' not in self.description):
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io, 'K_io': self.K_io,
                                 'K_io_f': self.K_io_f, 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.extended_flows_names, 'STR': self.STR_f.to_dict(),
                                 'impact_categories_IO': self.extended_impact_names_IW_exio,
                                 'impact_categories_eco': self.extended_impact_names_IW_eco}
            elif ('Environmental extensions were completed' in self.description and
                    'Impact World+ was used' not in self.description):
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io, 'K_io': self.K_io,
                                 'K_io_f': self.K_io_f, 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.extended_flows_names, 'STR': self.STR_f.to_dict(),
                                 'IMP': self.IMP.to_dict(), 'impact_categories_IO': self.extended_impact_names_CML}
            elif 'Regionalized flows/impacts available' in self.description:
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io, 'K_io': self.K_io,
                                 'K_io_f': self.K_io_f, 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.extended_flows_names, 'STR': self.STR_f.to_dict(),
                                 'impact_categories_IO': self.extended_impact_names_IW_exio,
                                 'impact_categories_eco': self.extended_impact_names_IW_eco, 'F_f_regio': self.F_f_regio,
                                 'C_io_regio': self.C_io_regio, 'F_io_regio': self.F_io_regio, 'C_f_regio': self.C_f_regio,
                                 'regionalized_flow_names_eco': self.regionalized_flow_names_eco,
                                 'regionalized_flow_names_exio': self.regionalized_flow_names_exio,
                                 'regionalized_impact_names_eco': self.regionalized_impact_names_eco,
                                 'regionalized_impact_names_exio': self.regionalized_impact_names_exio}
            else:
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io, 'K_io': self.K_io,
                                 'K_io_f': self.K_io_f, 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.flows_of_IO, 'impact_categories_IO': self.impact_methods_IO,
                                 'IMP': self.IMP.to_dict(), 'STR': self.STR_f.to_dict()}
        else:
            if ('Impact World+ was used' in self.description and
                    'Regionalized flows/impacts available' not in self.description):
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io,
                                 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.extended_flows_names, 'STR': self.STR_f.to_dict(),
                                 'impact_categories_IO': self.extended_impact_names_IW_exio,
                                 'impact_categories_eco': self.extended_impact_names_IW_eco}
            elif ('Environmental extensions were completed' in self.description and
                    'Impact World+ was used' not in self.description):
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io,
                                 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.extended_flows_names, 'STR': self.STR_f.to_dict(),
                                 'IMP': self.IMP.to_dict(), 'impact_categories_IO': self.extended_impact_names_CML}
            elif 'Regionalized flows/impacts available' in self.description:
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io,
                                 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.extended_flows_names, 'STR': self.STR_f.to_dict(),
                                 'impact_categories_IO': self.extended_impact_names_IW_exio,
                                 'impact_categories_eco': self.extended_impact_names_IW_eco, 'F_f_regio': self.F_f_regio,
                                 'C_io_regio': self.C_io_regio, 'F_io_regio': self.F_io_regio, 'C_f_regio': self.C_f_regio,
                                 'regionalized_flow_names_eco': self.regionalized_flow_names_eco,
                                 'regionalized_flow_names_exio': self.regionalized_flow_names_exio,
                                 'regionalized_impact_names_eco': self.regionalized_impact_names_eco,
                                 'regionalized_impact_names_exio': self.regionalized_impact_names_exio}
            else:
                hybrid_system = {'PRO_f': self.PRO_f.to_dict(), 'A_ff': self.A_ff, 'A_io': self.A_io, 'A_io_f': self.A_io_f,
                                 'F_f': self.F_f, 'F_io': self.F_io, 'C_f': self.C_f, 'C_io': self.C_io,
                                 'hybridized_processes': self.hybridized_processes, 'description': self.description,
                                 'sectors_of_IO': self.sectors_of_IO, 'regions_of_IO': self.regions_of_IO,
                                 'flows_of_IO': self.flows_of_IO, 'impact_categories_IO': self.impact_methods_IO,
                                 'IMP': self.IMP.to_dict(), 'STR': self.STR_f.to_dict()}

        if not os.path.exists(pkg_resources.resource_filename(
                __name__, '/Databases/' + self.lca_database_name_and_version + '_' + self.io_database_name_and_version)):
            os.makedirs(pkg_resources.resource_filename(
                __name__, '/Databases/' + self.lca_database_name_and_version + '_' + self.io_database_name_and_version))

        if not os.path.exists(pkg_resources.resource_filename(
                __name__, '/Databases/' + self.lca_database_name_and_version + '_' +
                          self.io_database_name_and_version + '/__init__.py')):
            os.makedirs(pkg_resources.resource_filename(
                __name__, '/Databases/' + self.lca_database_name_and_version + '_' +
                          self.io_database_name_and_version + '/__init__.py'))

        file = open(pkg_resources.resource_filename(
            __name__, '/Databases/' + self.lca_database_name_and_version + '_' +
                      self.io_database_name_and_version+'/description_system.txt'), 'w')
        file.write(str(self.description))
        file.close()

        if format == 'pickle':
            with gzip.open((pkg_resources.resource_filename(
                    __name__, '/Databases/' + self.lca_database_name_and_version + '_' +
                              self.io_database_name_and_version + '/hybrid_system.pickle')), 'wb') as f:
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

    def __init__(self, path_to_hybrid_system):

        try:
            with gzip.open(path_to_hybrid_system, 'rb') as f:
                self.hybrid_system = pd.read_pickle(f)
        except FileNotFoundError:
            print('The path entered does not lead to the required hybrid_system.pickle file.')
        try:
            self.description = eval(open(path_to_hybrid_system.split('/hybrid_system.pickle')[0]+
                                         '/description_system.txt').read())
        except FileNotFoundError:
            print('Please put the description_system.txt file in the same folder as your hybrid_system.pickle file')

        self.PRO_f = self.hybrid_system['PRO_f']
        self.A_ff = self.hybrid_system['A_ff']
        self.A_io = self.hybrid_system['A_io']
        self.A_io_f = self.hybrid_system['A_io_f']
        self.F_f = self.hybrid_system['F_f']
        self.F_io = self.hybrid_system['F_io']
        self.C_f = self.hybrid_system['C_f']
        self.C_io = self.hybrid_system['C_io']
        self.STR = self.hybrid_system['STR']
        self.hybridized_processes = self.hybrid_system['hybridized_processes']
        self.description = self.hybrid_system['description']
        self.flows_of_IO = self.hybrid_system['flows_of_IO']
        self.regions_of_IO = self.hybrid_system['regions_of_IO']
        self.sectors_of_IO = self.hybrid_system['sectors_of_IO']
        self.impact_categories_IO = self.hybrid_system['impact_categories_IO']
        try:
            self.K_io = self.hybrid_system['K_io']
            self.K_io_f = self.hybrid_system['K_io_f']
        except KeyError:
            pass
        if 'Impact World+ was used' not in self.description:
            self.IMP = self.hybrid_system['IMP']
            try :
                self.C = pd.concat([pd.DataFrame(self.C_f.todense(), index=self.IMP['method'], columns=self.STR['MATRIXID']),
                           pd.DataFrame(self.C_io.todense(), index=self.impact_categories_IO,
                                        columns=self.flows_of_IO)]).fillna(0)
            except KeyError:
                self.C = pd.concat([pd.DataFrame(self.C_f.todense(), index=self.IMP['Method'], columns=self.STR['MATRIXID']),
                           pd.DataFrame(self.C_io.todense(), index=self.impact_categories_IO,
                                        columns=self.flows_of_IO)]).fillna(0)
            self.C_index = self.C.index.tolist()
            self.C = back_to_sparse(self.C)
        else:
            self.impact_categories_eco = self.hybrid_system['impact_categories_eco']
            if 'Regionalized flows/impacts available' not in self.description:
                self.C = pd.DataFrame(self.C_f.todense(), self.impact_categories_eco,
                                      list(self.STR['MATRIXID'])).merge(
                    pd.DataFrame(self.C_io.todense(), self.impact_categories_IO, self.flows_of_IO),
                    left_on=pd.DataFrame(self.C_f.todense(), self.impact_categories_eco,
                                         list(self.STR['MATRIXID'])).index,
                    right_on=pd.DataFrame(self.C_io.todense(), self.impact_categories_IO,
                                          self.flows_of_IO).index,
                    how='outer')
                self.C_index = self.C['key_0'].tolist()
                self.C.drop(['key_0'], axis=1, inplace=True)
                self.C = self.C.fillna(0)
                self.C = back_to_sparse(self.C)
            else:
                self.C_io_regio = self.hybrid_system['C_io_regio']
                self.C_f_regio = self.hybrid_system['C_f_regio']
                self.F_io_regio = self.hybrid_system['F_io_regio']
                self.F_f_regio = self.hybrid_system['F_f_regio']
                self.regionalized_impact_names_exio = self.hybrid_system['regionalized_impact_names_exio']
                self.regionalized_impact_names_eco = self.hybrid_system['regionalized_impact_names_eco']
                self.regionalized_flow_names_exio = self.hybrid_system['regionalized_flow_names_exio']
                self.regionalized_flow_names_eco = self.hybrid_system['regionalized_flow_names_eco']

                self.C = pd.concat([pd.DataFrame(self.C_f.todense(), self.impact_categories_eco, self.STR['MATRIXID']),
                                    pd.DataFrame(self.C_io.todense(), self.impact_categories_IO, self.flows_of_IO)],
                                   sort=False).fillna(0)
                self.C = self.C.groupby(self.C.index).sum()
                self.non_regio_impacts_index = self.C.index.tolist()
                self.C = back_to_sparse(self.C)
                self.C_regio = pd.concat([pd.DataFrame(self.C_f_regio.todense(), index=self.regionalized_flow_names_eco,
                                                       columns=self.regionalized_impact_names_eco),
                                          pd.DataFrame(self.C_io_regio.todense(), index=self.regionalized_flow_names_exio,
                                                       columns=self.regionalized_impact_names_exio)],
                                         sort=False).fillna(0).T
                self.C_regio = self.C_regio.groupby(self.C_regio.index).sum()
                self.regio_impacts_index = self.C_regio.index.tolist()
                self.C_regio = back_to_sparse(self.C_regio)

        del self.hybrid_system

        self.D = pd.DataFrame()
        self.Lp = pd.DataFrame()
        self.Lio = pd.DataFrame()
        self.original = pd.DataFrame()
        self.total_emissions_not_regio = pd.DataFrame()
        self.total_emissions_regio = pd.DataFrame()
        self.total_emissions = pd.DataFrame()

        self.GWP100_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/GWP.txt').decode('utf-8'))
        self.Acidification_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/Acid.txt').decode('utf-8'))
        self.Eutrophication_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/Eutro.txt').decode('utf-8'))
        self.HTox_CML2001 = ast.literal_eval(
            pkg_resources.resource_string(__name__, '/Data/Characterization_matching/HTox.txt').decode('utf-8'))

    def calc_lifecycle(self):
        """ Simply calculates lifecycle production, emissions, or impacts

        Returns
        -------
            * lifecycle impact in self.D

        """

        self.Lp = scipy.sparse.csr_matrix(
            np.linalg.solve(np.eye(self.A_ff.shape[0]) - self.A_ff.todense(), np.eye(self.A_ff.shape[0])))

        if 'Capitals were endogenized' in self.description:
            self.Lio = scipy.sparse.csr_matrix(
                np.linalg.solve(np.eye(self.A_io.shape[0]) - (self.A_io+self.K_io).todense(), np.eye(self.A_io.shape[0])))
            if 'Regionalized flows/impacts available' in self.description:
                self.total_emissions_not_regio = scipy.sparse.vstack(
                    [self.F_f.dot(self.Lp), self.F_io.dot(self.Lio).dot(self.A_io_f + self.K_io_f).dot(self.Lp)])
                self.total_emissions_regio = scipy.sparse.vstack([self.F_f_regio.dot(self.Lp),
                                                                  self.F_io_regio.dot(self.Lio).dot(
                                                                      self.A_io_f+self.K_io_f).dot(self.Lp)])
                self.D = pd.DataFrame(scipy.sparse.vstack([self.C.dot(self.total_emissions_not_regio),
                                                           self.C_regio.dot(self.total_emissions_regio)]).todense(),
                                      index=self.non_regio_impacts_index+self.regio_impacts_index,
                                      columns=self.PRO_f['activityId']).loc[:, self.hybridized_processes]
            else:
                self.total_emissions = scipy.sparse.vstack([self.F_f.dot(self.Lp), self.F_io.dot(self.Lio).dot(
                    self.A_io_f+self.K_io_f).dot(self.Lp)])
                self.D = pd.DataFrame(
                    (self.C.dot(self.total_emissions)).todense(),
                    index=self.C_index, columns=self.PRO_f['activityName']).loc[:, self.hybridized_processes]

            print('Calculations done! Results are contained in self.D')

        else:
            self.Lio = scipy.sparse.csr_matrix(
                np.linalg.solve(np.eye(self.A_io.shape[0]) - self.A_io.todense(), np.eye(self.A_io.shape[0])))
            if 'Regionalized flows/impacts available' in self.description:
                self.total_emissions_not_regio = scipy.sparse.vstack(
                    [self.F_f.dot(self.Lp), self.F_io.dot(self.Lio).dot(self.A_io_f).dot(self.Lp)])
                self.total_emissions_regio = scipy.sparse.vstack([self.F_f_regio.dot(self.Lp),
                                                                  self.F_io_regio.dot(self.Lio).dot(
                                                                      self.A_io_f).dot(self.Lp)])
                self.D = pd.DataFrame(scipy.sparse.vstack([self.C.dot(self.total_emissions_not_regio),
                                                           self.C_regio.dot(self.total_emissions_regio)]).todense(),
                                      index=self.non_regio_impacts_index+self.regio_impacts_index,
                                      columns=self.PRO_f['activityId']).loc[:, self.hybridized_processes]
            else:
                self.total_emissions = scipy.sparse.vstack([self.F_f.dot(self.Lp), self.F_io.dot(self.Lio).dot(
                    self.A_io_f).dot(self.Lp)])
                self.D = pd.DataFrame(
                    (self.C.dot(self.total_emissions)).todense(),
                    index=self.C_index, columns=self.PRO_f['activityName']).loc[:, self.hybridized_processes]
            print('Calculations done! Results are contained in self.D')

    def contributions(self, process):

        try:
            if self.Lp.empty:
                self.Lp = scipy.sparse.csr_matrix(
                    np.linalg.solve(np.eye(self.A_ff.shape[0]) - self.A_ff.todense(), np.eye(self.A_ff.shape[0])))
            if self.Lio.empty:
                if 'Capitals were endogenized' in self.description:
                    self.Lio = scipy.sparse.csr_matrix(
                        np.linalg.solve(np.eye(self.A_io.shape[0]) - (self.A_io + self.K_io).todense(),
                                        np.eye(self.A_io.shape[0])))
                else:
                    self.Lio = scipy.sparse.csr_matrix(
                        np.linalg.solve(np.eye(self.A_io.shape[0]) - self.A_io.todense(), np.eye(self.A_io.shape[0])))
        except AttributeError:
            pass

        position_process = [i for i, x in enumerate(list(self.PRO_f['activityId'])) if x == process][0]

        Xp = self.Lp.dot(scipy.sparse.csr_matrix(np.diagflat(self.A_ff.todense()[:, position_process])))

        if 'Capitals were endogenized' in self.description:
            Xiop = self.Lio.dot(self.A_io_f+self.K_io_f).dot(Xp)

            Xio = self.Lio.dot(scipy.sparse.csr_matrix(np.diagflat((self.A_io_f+self.K_io_f).todense()[:, position_process])))
        else:
            Xiop = self.Lio.dot(self.A_io_f).dot(Xp)

            Xio = self.Lio.dot(scipy.sparse.csr_matrix(np.diagflat(self.A_io_f.todense()[:, position_process])))

        X = scipy.sparse.hstack([scipy.sparse.vstack([Xp, Xiop]), scipy.sparse.vstack(
            [scipy.sparse.csr_matrix(np.zeros((Xp.shape[0], Xio.shape[0]))), Xio])])

        F = scipy.sparse.hstack([scipy.sparse.vstack(
            [self.F_f, scipy.sparse.csr_matrix(np.zeros((self.F_io.shape[0], self.A_ff.shape[0])))]),
                             scipy.sparse.vstack(
                                 [scipy.sparse.csr_matrix(np.zeros((self.F_f.shape[0], self.A_io.shape[0]))),
                                  self.F_io])])

        full_index = list(self.PRO_f['activityId'])+pd.MultiIndex.from_product([
            self.regions_of_IO, self.sectors_of_IO]).tolist()

        if 'Regionalized flows/impacts available' in self.description:
            F_regio = scipy.sparse.hstack([scipy.sparse.vstack(
                [self.F_f_regio, scipy.sparse.csr_matrix(np.zeros((self.F_io_regio.shape[0], self.A_ff.shape[0])))]),
                scipy.sparse.vstack(
                    [scipy.sparse.csr_matrix(np.zeros((self.F_f_regio.shape[0], self.A_io.shape[0]))),
                     self.F_io_regio])])

            return (pd.DataFrame(X.todense(), index=full_index, columns=full_index),
                    pd.DataFrame(F.todense(), index=list(self.STR['MATRIXID'])+self.flows_of_IO, columns=full_index),
                    pd.DataFrame(self.C.todense(), index=self.non_regio_impacts_index,
                                 columns=list(self.STR['MATRIXID']) + self.flows_of_IO),
                    pd.DataFrame(F_regio.todense(),
                                 index=self.regionalized_flow_names_eco+self.regionalized_flow_names_exio,
                                 columns=full_index),
                    pd.DataFrame(self.C_regio.todense(),index=self.regio_impacts_index,
                                 columns=self.regionalized_flow_names_eco+self.regionalized_flow_names_exio))

        else:
            return (pd.DataFrame(X.todense(), index=full_index, columns=full_index),
                    pd.DataFrame(F.todense(), index=list(self.STR['MATRIXID'])+self.flows_of_IO, columns=full_index),
                    pd.DataFrame(self.C.todense(), index=self.C_index,
                                 columns=list(self.STR['MATRIXID'])+self.flows_of_IO))

    def get_available_impact_methods(self, impact_category):

        if impact_category == 'GWP':
            return self.GWP100_CML2001
        elif impact_category == 'Acidification':
            return self.Acidification_CML2001
        elif impact_category == 'Eutrophication':
            return self.Eutrophication_CML2001
        elif impact_category == 'Human toxicity':
            return self.HTox_CML2001
        else:
            print('Possible impact categories: "GWP", "Acidification", "Eutrophication", "Human toxicity"')

    def navigate_through_PRO_f(self, product=False, geography=False, activity=False):

        PRO_f = pd.DataFrame.from_dict(self.PRO_f)

        if product and geography and activity:
            return PRO_f.loc[[i for i in PRO_f.index if product in PRO_f.productName[i] and
                                   geography in PRO_f.geography[i] and
                                   activity in PRO_f.activityName[i]]]
        elif product and geography:
            return PRO_f.loc[[i for i in PRO_f.index if product in PRO_f.productName[i] and
                                   geography in PRO_f.geography[i]]]
        elif activity and geography:
            return PRO_f.loc[[i for i in PRO_f.index if activity in PRO_f.activityName[i] and
                                   geography in PRO_f.geography[i]]]
        elif product:
            return PRO_f.loc[[i for i in PRO_f.index if product in PRO_f.productName[i]]]
        elif activity:
            return PRO_f.loc[[i for i in PRO_f.index if activity in PRO_f.activityName[i]]]
        else:
            print('Enter at least a product or an activity')

    def export_to_brightway2(self, bw2_project_name, created_database_name, path_to_ecoinvent_ecospold_datasets,
                             aggregated=False):
        """
        Exporting the hybrid database into brightway2 in a disaggregated format.

        :args:  bw2_project_name: (string) the name of the brightway2 project in which the hybrid database will be
                                           imported. Note that the method creates and deletes an 'hybrid-ecoinvent'
                                           database. If you use an already-existing project and have a database named
                                           'hybrid-ecoinvent', it will be overwritten and deleted.
               path_to_ecoinvent_ecospold_datasets: (string) path to the ecospold files of the ecoinvent version used
               aggregated: (boolean) if True Contribution analyses will show the contribution aggregated IO complements
                                             comparatively to ecoinvent inputs. Takes 1h30 to 2h to import.
                                     if False Contribution analyses will show the contribution of each IO complement
                                              comparatively to ecoinvent inputs. Takes 4 to 5h to import.

        :return: nothing, but it creates two databases ('exiobase' and 'hybrid_ecoinvent') into the selected
                 brightway2 project
        """

        if 'Regionalized flows/impacts available' in self.description:
            print('Cannot export regionalized version to brightway2.')
            return

        # creating a new bw2 project
        projects.set_current(bw2_project_name)

        # importing elementary flows and default LCIA methods
        bw2setup()

        # loading the IMPACT World+ LCIA method into the project
        BW2Package.import_file(pkg_resources.resource_stream(__name__, '/Data/for_export_in_brightway2/'
                                                                       'laurepatou-IMPACT-World-in-Brightway-d949b33/'
                                                                       'Brightway_IW_damage_1_46_and_midpoint_1_28.bw2package'))

        # importing hybrid-ecoinvent into the project
        eco_importer = SingleOutputEcospold2Importer(path_to_ecoinvent_ecospold_datasets, 'hybrid-ecoinvent')
        eco_importer.apply_strategies()
        eco_importer.write_database()

        if not aggregated:
            # importing exiobase into the project
            self.import_exiobase_into_brightway2()

        # importing hybrid ecoinvent
        self.import_hybrid_ecoinvent_into_brightway2(created_database_name, aggregated)

    def import_exiobase_into_brightway2(self):

        print('Starting import of exiobase')
        print('Aggregating results')
        # calculating aggregated results for exiobase
        self.Lio = scipy.sparse.csr_matrix(
            np.linalg.solve(np.eye(self.A_io.shape[0]) - self.A_io.todense(), np.eye(self.A_io.shape[0])))
        aggregated_results = pd.DataFrame((self.C_io.dot(self.F_io).dot(self.Lio)).todense(), self.impact_categories_IO,
                                          pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO]))
        print('Formating data')
        # create hash codes for exiobase sectors
        give_hash_to_exio = dict.fromkeys(aggregated_results.columns)
        give_hash_to_exio = {k: hashlib.md5(str(k).encode()).hexdigest() for k, v in give_hash_to_exio.items()}
        aggregated_results.columns = [give_hash_to_exio[i] for i in aggregated_results.columns]

        # link IMPACT World+ impact categories typology to the methods abbreviations used in brightway2
        IW_pylcaio_to_bw2 = dict.fromkeys(self.impact_categories_IO)
        for IW_category in IW_pylcaio_to_bw2:
            if 'PDF' in IW_category or 'DALY' in IW_category:
                IW_pylcaio_to_bw2[IW_category] = \
                methods.get(('IMPACTWorld+ (Default_Recommended_Damage 1.46)', IW_category.split(' (')[0]))[
                    'abbreviation']
            else:
                IW_pylcaio_to_bw2[IW_category] = \
                methods.get(('IMPACTWorld+ (Default_Recommended_Midpoint 1.28)', IW_category.split(' (')[0]))[
                    'abbreviation']
        aggregated_results.index = list(IW_pylcaio_to_bw2.values())

        # create dummy exchanges and characterization factors to enable the storing of LCIA scores into exchanges
        add_all_unit_score_exchanges_and_cfs()

        # formatting exiobase data
        bw2_dict = dict.fromkeys(list(zip(['exiobase'] * len(aggregated_results.columns), aggregated_results.columns)),
                                 {
                                     'name': '',
                                     'unit': '',
                                     'location': '',
                                     'exchanges': [],
                                 })

        for sector_hash in aggregated_results.columns:
            bw2_dict[('exiobase', sector_hash)] = {
                'name': list({k: v for k, v in give_hash_to_exio.items() if v == sector_hash}.keys())[0][1],
                'unit': 'euro',
                'location': list({k: v for k, v in give_hash_to_exio.items() if v == sector_hash}.keys())[0][0],
                'exchanges': aggregated_results.loc[:, sector_hash][
                    aggregated_results.loc[:, sector_hash] != 0].to_dict()}
        for sector_hash in aggregated_results.columns:
            bw2_dict[('exiobase', sector_hash)]['exchanges'] = list(
                {k: {'input': ('biosphere3', k), 'amount': v, 'type': 'biosphere'} for k, v in
                 bw2_dict[('exiobase', sector_hash)]['exchanges'].items()}.values())
        print('Writing exiobase database')
        # importing the aggregated exiobase into the project
        Database('exiobase').write(bw2_dict)
        print('Created exiobase database')

    def import_hybrid_ecoinvent_into_brightway2(self, created_database_name, aggregated):
        print('Starting import of hybrid data')
        db = Database('hybrid-ecoinvent')
        ecoinvent_metadata = pd.DataFrame.from_dict(self.PRO_f)

        # need to match uuids used in pylcaio to hash codes used by brightway2
        match_hash_to_uuid = dict.fromkeys([act.as_dict()['code'] for act in db])
        for act in db:
            match_hash_to_uuid[act.as_dict()['code']] = act.as_dict()['activity'] + '_' + act.as_dict()['flow']

        if aggregated:
            print('Aggregating data')
            aggregated_results = self.aggregate_hybrid_ecoinvent()
            print('Formating hybrid data')
            # match LCIA methods
            IW_pylcaio_to_bw2 = dict.fromkeys(aggregated_results.index)
            for IW_category in IW_pylcaio_to_bw2:
                if 'PDF' in IW_category or 'DALY' in IW_category:
                    IW_pylcaio_to_bw2[IW_category] = \
                    methods.get(('IMPACTWorld+ (Default_Recommended_Damage 1.46)', IW_category.split(' (')[0]))[
                        'abbreviation']
                else:
                    IW_pylcaio_to_bw2[IW_category] = \
                    methods.get(('IMPACTWorld+ (Default_Recommended_Midpoint 1.28)', IW_category.split(' (')[0]))[
                        'abbreviation']
            aggregated_results.index = list(IW_pylcaio_to_bw2.values())

            match_uuid_to_hash = {v: k for k, v in match_hash_to_uuid.items()}
            aggregated_results.columns = [match_uuid_to_hash[i] for i in aggregated_results.columns]

            add_all_unit_score_exchanges_and_cfs()

            bw2_dict = dict.fromkeys(list(zip([created_database_name] * len(aggregated_results.columns),
                                              aggregated_results.columns)), {
                'exchanges': [],
            })

            for process_hash in aggregated_results.columns:
                bw2_dict[(created_database_name, process_hash)] = {
                    'exchanges': aggregated_results.loc[:, process_hash][aggregated_results.loc[:, process_hash] != 0].to_dict()}
            for process_hash in aggregated_results.columns:
                bw2_dict[(created_database_name, process_hash)]['exchanges'] = list(
                    {k: {'input': ('biosphere3', k), 'amount': v, 'type': 'biosphere'} for k, v in
                     bw2_dict[(created_database_name, process_hash)]['exchanges'].items()}.values())
            print('Matching ecoinvent data with hybrid data')
            # to link ecoinvent exchanges to exiobase exchanges, we need the original ecoinvent to also be in a dict format
            bw2_eco_dict = db.load()

            # add IO complements to each ecoinvent process
            for process in bw2_eco_dict.keys():
                bw2_eco_dict[process]['exchanges'].extend(bw2_dict[(created_database_name, process[1])]['exchanges'])

            # exchanges will be associated to the hybrid_ecoinvent_agg database of bw2 and not hybrid-ecoinvent
            bw2_hybrid_dict = {(created_database_name, k[1]): v for k, v in bw2_eco_dict.items()}
            bw2_hybrid_dict = {
                outer_k: {
                    inner_k: [
                        {k: (created_database_name, v[1]) if ((k == 'input' or k == 'output')
                                                           and v[0] == 'hybrid-ecoinvent')
                         else v
                         for k, v in i.items()}
                        for i in bw2_hybrid_dict[outer_k][inner_k]
                    ]
                    if inner_k == 'exchanges' else inner_v
                    for inner_k, inner_v in outer_v.items()
                }
                for outer_k, outer_v in bw2_hybrid_dict.items()
            }
            print('Deleting hybrid-ecoinvent database')
            del databases['hybrid-ecoinvent']
            print('Writing {} database'.format(created_database_name))
            Database(created_database_name).write(bw2_hybrid_dict)
            print('{} successfully imported'.format(created_database_name))

        else:
            A_io_f = pd.DataFrame(self.A_io_f.todense(),
                                  pd.MultiIndex.from_product([self.regions_of_IO, self.sectors_of_IO]),
                                  ecoinvent_metadata.index)
            print('Formating hybrid data')
            # create hash codes for exiobase sectors
            give_hash_to_exio = dict.fromkeys(A_io_f.index)
            give_hash_to_exio = {k: hashlib.md5(str(k).encode()).hexdigest() for k, v in give_hash_to_exio.items()}
            # set those hash codes as index
            A_io_f.index = [give_hash_to_exio[i] for i in list(A_io_f.index)]

            # format the results from pylcaio (the A_io_f matrix) to be readable by brightway2
            bw2_dict = dict.fromkeys(list(zip([created_database_name] * len(ecoinvent_metadata.index), match_hash_to_uuid)), {
                'name': '',
                'unit': '',
                'location': '',
                'exchanges': [],
            })
            for ecoinvent_uuid in ecoinvent_metadata.index:
                hash_code = list({k: v for k, v in match_hash_to_uuid.items() if v == ecoinvent_uuid}.keys())[0]
                bw2_dict[(created_database_name, hash_code)] = {'name': ecoinvent_metadata.loc[ecoinvent_uuid, 'activityName'],
                                                             'unit': ecoinvent_metadata.loc[ecoinvent_uuid, 'unitName'],
                                                             'location': ecoinvent_metadata.loc[
                                                                 ecoinvent_uuid, 'geography'],
                                                             'exchanges': A_io_f.loc[:, ecoinvent_uuid][
                                                                 A_io_f.loc[:, ecoinvent_uuid] != 0].to_dict()}
            for ecoinvent_uuid in ecoinvent_metadata.index:
                hash_code = list({k: v for k, v in match_hash_to_uuid.items() if v == ecoinvent_uuid}.keys())[0]
                bw2_dict[(created_database_name, hash_code)]['exchanges'] = list(
                    {k: {'input': ('exiobase', k), 'amount': v, 'type': 'technosphere'} for k, v in
                     bw2_dict[(created_database_name, hash_code)]['exchanges'].items()}.values())
            print('Matching ecoinvent data with hybrid data')
            # to link ecoinvent exchanges to exiobase exchanges, we need the original ecoinvent to also be in a dict format
            bw2_eco_dict = Database('hybrid-ecoinvent').load()

            # add IO complements to each ecoinvent process
            for process in bw2_eco_dict.keys():
                bw2_eco_dict[process]['exchanges'].extend(bw2_dict[(created_database_name, process[1])]['exchanges'])

            # exchanges will be associated to the hybrid_ecoinvent database of bw2 and not hybrid-ecoinvent
            bw2_hybrid_dict = {(created_database_name, k[1]): v for k, v in bw2_eco_dict.items()}

            # we also need to replace the database name in the inputs. It's a bit trickier:
            bw2_hybrid_dict = {
                outer_k: {
                    inner_k: [
                        {k: (created_database_name, v[1]) if ((k == 'input' or k == 'output')
                                                           and v[0] == 'hybrid-ecoinvent')
                         else v
                         for k, v in i.items()}
                        for i in bw2_hybrid_dict[outer_k][inner_k]
                    ]
                    if inner_k == 'exchanges' else inner_v
                    for inner_k, inner_v in outer_v.items()
                }
                for outer_k, outer_v in bw2_hybrid_dict.items()
            }
            print('Deleting hybrid-ecoinvent database')
            # we remove the 'hybrid-ecoinvent' database. Its taking unnecessary space.
            del databases['hybrid-ecoinvent']
            print('Writing {} database'.format(created_database_name))
            # write the formatted dictionary in the brightway2 database
            Database(created_database_name).write(bw2_hybrid_dict)
            print('{} successfully imported'.format(created_database_name))

    def aggregate_hybrid_ecoinvent(self):

        if 'Capitals were endogenized' in self.description:
            self.Lio = scipy.sparse.csr_matrix(
                np.linalg.solve(np.eye(self.A_io.shape[0]) - (self.A_io+self.K_io).todense(), np.eye(self.A_io.shape[0])))

        else:
            self.Lio = scipy.sparse.csr_matrix(
                np.linalg.solve(np.eye(self.A_io.shape[0]) - self.A_io.todense(), np.eye(self.A_io.shape[0])))

        impact_hybrid = pd.DataFrame(((self.C_io.dot(self.F_io).dot(self.Lio)).dot(self.A_io_f)).todense(),
                                     index=self.impact_categories_IO, columns=self.PRO_f['activityName'])
        return impact_hybrid


def get_inflation(reference_year):
    """ Returns the inflation rate between the year 2005 (base year for ecoinvent prices) and the reference year of
    the used IO database, from https://www.inflationtool.com/euro/2005-to-present-value"""

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
    elif reference_year == 2012:
        inflation = 1.16
    elif reference_year == 2013:
        inflation = 1.18
    elif reference_year == 2014:
        inflation = 1.19
    elif reference_year == 2015:
        inflation = 1.19
    elif reference_year == 2016:
        inflation = 1.19
    elif reference_year == 2017:
        inflation = 1.21
    elif reference_year == 2018:
        inflation = 1.22
    elif reference_year == 2019:
        inflation = 1.24
    elif reference_year == 2020:
        inflation = 1.26
    elif reference_year == 2021:
        inflation = 1.25
    # no data available for 2022, same data as 2021 by default
    elif reference_year == 2022:
        inflation = 1.25
    else:
        inflation = 1

    return inflation


def sum_elements_list(liste):
    concatenated_list = []
    for i in range(0, len(liste)):
        if isinstance(liste[i], list):
            concatenated_list += liste[i]
        else:
            concatenated_list += [liste[i]]
    return concatenated_list


def LCA_convention_to_IO(dataframe):
    """ Changes the convention of a technology matrix from LCA to IO """
    # no ones on the diagonal
    dataframe[dataframe == 1] = 0
    # only positive values
    dataframe = dataframe.abs()
    return dataframe


def back_to_sparse(dataframe):
    return scipy.sparse.csr_matrix(dataframe)


def completing_extensions(OG_extensions, new_extensions):
    """Function to modify the names of the extensions of the original exiobase to match the ones resulting from the
    matching with USEEIO. Also concatenates both original and new extensions at the end, resulting in the extended
    exiensions."""
    # just remove Energy Carrier Net flows and not defined waste flows
    OG_extensions.drop([i for i in OG_extensions.index if 'Energy Carrier Net' in i],inplace=True)
    OG_extensions.drop(['Emissions nec - waste - undef'],inplace=True)
    # simple sum of all values for old names pollutants, e.g., 'CH4 - non combustion - Cement production - air' => 'CH4 - air'
    easy_match = ['CH4','N2O','SOx','NH3','HCB','NMVOC','PM10','PM2.5','TSP','Cd','Hg','Pb','Zn','PAH']
    for pollutant in easy_match:
        list_old_names = [i for i in OG_extensions.index if pollutant in i]
        OG_extensions.loc[pollutant+' - air'] = OG_extensions.loc[list_old_names].sum()
        OG_extensions.drop(list_old_names, inplace=True)
    # hardcoded stuff
    OG_extensions.loc['Benzo(a)pyrene - air'] = OG_extensions.loc[['Benzo(a)pyrene - combustion - air','B(a)P - non combustion - Primary aluminium production - air','B(a)P - non combustion - Production of coke oven coke - air','B(a)P - non combustion - Production of gascoke - air']].sum()
    OG_extensions.drop(['Benzo(a)pyrene - combustion - air','B(a)P - non combustion - Primary aluminium production - air','B(a)P - non combustion - Production of coke oven coke - air','B(a)P - non combustion - Production of gascoke - air'],inplace=True)
    OG_extensions.loc['Benzo(b)fluoranthene - air'] = OG_extensions.loc[['Benzo(b)fluoranthene - combustion - air','B(b)F - non combustion - Primary aluminium production - air','B(b)F - non combustion - Production of coke oven coke - air','B(b)F - non combustion - Production of gascoke - air']].sum()
    OG_extensions.drop(['Benzo(b)fluoranthene - combustion - air','B(b)F - non combustion - Primary aluminium production - air','B(b)F - non combustion - Production of coke oven coke - air','B(b)F - non combustion - Production of gascoke - air'],inplace=True)
    OG_extensions.loc['Benzo(k)fluoranthene - air'] = OG_extensions.loc[['Benzo(k)fluoranthene - combustion - air','B(k)F - non combustion - Primary aluminium production - air','B(k)F - non combustion - Production of coke oven coke - air','B(k)F - non combustion - Production of gascoke - air']].sum()
    OG_extensions.drop(['Benzo(k)fluoranthene - combustion - air','B(k)F - non combustion - Primary aluminium production - air','B(k)F - non combustion - Production of coke oven coke - air','B(k)F - non combustion - Production of gascoke - air'],inplace=True)
    OG_extensions.loc['CO2 - biogenic - air'] = OG_extensions.loc['CO2 - waste - biogenic - air']
    OG_extensions.drop(['CO2 - waste - biogenic - air'],inplace=True)
    OG_extensions.loc['Pxx - soil'] = OG_extensions.loc[['Pxx - agriculture - soil','P - agriculture - soil']].sum()
    OG_extensions.drop(['Pxx - agriculture - soil','P - agriculture - soil'],inplace=True)
    old_CO2_flows = [i for i in OG_extensions.index if 'CO2' in i and 'biogenic' not in i]
    OG_extensions.loc['CO2 - air'] = OG_extensions.loc[old_CO2_flows].sum()
    OG_extensions.drop(old_CO2_flows,inplace=True)
    old_CO_flows = [i for i in OG_extensions.index if 'CO' in i and '2' not in i]
    OG_extensions.loc['CO - air'] = OG_extensions.loc[old_CO_flows].sum()
    OG_extensions.drop(old_CO_flows,inplace=True)
    old_NOx_flows = [i for i in OG_extensions.index if 'NOx' in i or 'NOX' in i]
    OG_extensions.loc['NOx - air'] = OG_extensions.loc[old_NOx_flows].sum()
    OG_extensions.drop(old_NOx_flows,inplace=True)
    old_indeno_flows = [i for i in OG_extensions.index if 'Indeno' in i]
    OG_extensions.loc['Indeno(1,2,3-cd)pyrene - air'] = OG_extensions.loc[old_indeno_flows].sum()
    OG_extensions.drop(old_indeno_flows,inplace=True)
    old_PCB_flows = [i for i in OG_extensions.index if 'PCB' in i]
    OG_extensions.loc['PCBs - air'] = OG_extensions.loc[old_PCB_flows].sum()
    OG_extensions.drop(old_PCB_flows,inplace=True)
    old_PCDD_flows = [i for i in OG_extensions.index if 'PCDD' in i]
    OG_extensions.loc['PCDD_F - air'] = OG_extensions.loc[old_PCDD_flows].sum()
    OG_extensions.drop(old_PCDD_flows,inplace=True)
    old_As_flows = [i for i in OG_extensions.index if 'As -' in i]
    OG_extensions.loc['As - air'] = OG_extensions.loc[old_As_flows].sum()
    OG_extensions.drop(old_As_flows,inplace=True)
    old_Ni_flows = [i for i in OG_extensions.index if 'Ni -' in i]
    OG_extensions.loc['Ni - air'] = OG_extensions.loc[old_Ni_flows].sum()
    OG_extensions.drop(old_Ni_flows,inplace=True)
    old_Cr_flows = [i for i in OG_extensions.index if 'Cr -' in i]
    OG_extensions.loc['Cr - air'] = OG_extensions.loc[old_Cr_flows].sum()
    OG_extensions.drop(old_Cr_flows,inplace=True)
    old_Cu_flows = [i for i in OG_extensions.index if 'Cu -' in i]
    OG_extensions.loc['Cu - air'] = OG_extensions.loc[old_Cu_flows].sum()
    OG_extensions.drop(old_Cu_flows,inplace=True)
    old_Se_flows = [i for i in OG_extensions.index if 'Se -' in i]
    OG_extensions.loc['Se - air'] = OG_extensions.loc[old_Se_flows].sum()
    OG_extensions.drop(old_Se_flows,inplace=True)
    old_N_flows = [i for i in OG_extensions.index if 'N -' in i]
    OG_extensions.loc['N - water'] = OG_extensions.loc[old_N_flows].sum()
    OG_extensions.drop(old_N_flows,inplace=True)
    old_P_flows = [i for i in OG_extensions.index if 'P -' in i and 'water' in i]
    OG_extensions.loc['P - water'] = OG_extensions.loc[old_P_flows].sum()
    OG_extensions.drop(old_P_flows,inplace=True)
    # after all this hardwork, concatenate with extensions
    extended_extensions = pd.concat([OG_extensions,new_extensions])
    return extended_extensions
