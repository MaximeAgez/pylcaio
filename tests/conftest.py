import numpy as np
import pandas as pd
import pytest
import uuid

@pytest.fixture()
def random_LCA_system():
    lca_database_name_and_version = 'ecoinvent3.5'
    activityId1 = str(uuid.uuid4())
    activityId2 = str(uuid.uuid4())
    activityId3 = str(uuid.uuid4())
    productId1 = str(uuid.uuid4())
    productId2 = str(uuid.uuid4())
    productId3 = str(uuid.uuid4())
    productId4 = str(uuid.uuid4())
    productId5 = str(uuid.uuid4())
    PRO_f = pd.DataFrame([['product1', productId1, 'activity1', activityId1, 'CA-QC', 'CA', 1.1],
                          ['product2', productId2, 'activity2', activityId2, 'RER', 'RER', 1.2],
                          ['product3', productId3, 'activity2', activityId2, 'RoW', 'RoW', 1.3],
                          ['product4', productId4, 'activity3', activityId3, 'CH', 'CH', 1.4],
                          ['product5', productId5, 'activity1', activityId1, 'RoW', 'RoW', 1.5]],
                         index=[str(productId1) + '_' + str(activityId1), str(productId2) + '_' + str(activityId2),
                                str(productId3) + '_' + str(activityId3),
                                str(productId4) + '_' + str(activityId3), str(productId5) + '_' + str(activityId1)],
                         columns=['productName', 'productId', 'activityName', 'activityId', 'geography', 'io_geography',
                                  'price'])
    A_ff = pd.DataFrame(np.random.random_sample((5, 5)),
                        index=[str(productId1) + '_' + str(activityId1), str(productId2) + '_' + str(activityId2),
                               str(productId3) + '_' + str(activityId3),
                               str(productId4) + '_' + str(activityId3), str(productId5) + '_' + str(activityId1)],
                        columns=[str(productId1) + '_' + str(activityId1), str(productId2) + '_' + str(activityId2),
                                 str(productId3) + '_' + str(activityId3),
                                 str(productId4) + '_' + str(activityId3),
                                 str(productId5) + '_' + str(activityId1)])

    F_f = pd.DataFrame(np.random.random_sample((3, 5)), index=['pollutant1', 'pollutant2', 'pollutant3'],
                       columns=[str(productId1) + '_' + str(activityId1), str(productId2) + '_' + str(activityId2),
                                str(productId3) + '_' + str(activityId3),
                                str(productId4) + '_' + str(activityId3), str(productId5) + '_' + str(activityId1)])
    C_f = pd.DataFrame(np.random.random_sample((3, 3)), index=['impact1', 'impact2', 'impact3'],
                       columns=['pollutant1', 'pollutant2', 'pollutant3'])
    STR_f = pd.DataFrame([['pollutant_name1', 'pollutant_unit1', 'pollutant_comp1', 'cas1'],
                          ['pollutant_name2', 'pollutant_unit2', 'pollutant_comp2', 'cas2'],
                          ['pollutant_name3', 'pollutant_unit3', 'pollutant_comp3', 'cas3']],
                         index=['pollutant1', 'pollutant2', 'pollutant3'],
                         columns=['pollutant_name', 'pollutant_unit', 'pollutant_comp', 'cas'])
    y_f = pd.DataFrame(np.eye(len(A_ff)), A_ff.index, A_ff.columns)

    dict_dataframes = {'PRO_f': PRO_f, 'A_ff': A_ff, 'F_f': F_f, 'C_f': C_f, 'STR_f': STR_f, 'y_f': y_f,
                       'lca_database_name_and_version': lca_database_name_and_version}
    return dict_dataframes

@pytest.fixture()
def random_IO_system():
    io_database_name_and_version = 'exiobase3'
    arrays = [['CA', 'FR', 'CH', 'ZA'], ['sector1', 'sector2', 'sector3', 'sector4', 'sector5']]
    index = pd.MultiIndex.from_product(arrays, names=['region', 'sector'])
    Z = pd.DataFrame(np.random.random_sample((20, 20)), index=index, columns=index)
    A_io = Z.mul(1 / (1.1 * Z.sum(axis=0)))
    F_io = pd.DataFrame(np.random.random_sample((3, 20)), index=['pollutant1.1', 'pollutant2.1', 'pollutant3.1'],
                        columns=index)
    C_io = pd.DataFrame(np.random.random_sample((3, 3)), index=['impact1.1', 'impact2.1', 'impact3.1'],
                        columns=['pollutant1.1', 'pollutant2.1', 'pollutant3.1'])
    STR_io = pd.DataFrame([['pollutant_name1.1', 'pollutant_unit1.1', 'pollutant_comp1.1', 'cas1.1'],
                           ['pollutant_name2.1', 'pollutant_unit2.1', 'pollutant_comp2.1', 'cas2.1'],
                           ['pollutant_name3.1', 'pollutant_unit3.1', 'pollutant_comp3.1', 'cas3.1']],
                          index=['pollutant1.1', 'pollutant2.1', 'pollutant3.1'],
                          columns=['pollutant_name', 'pollutant_unit', 'pollutant_comp', 'cas'])
    arraysFD = [['CA', 'FR', 'CH', 'ZA'], ['household', 'government']]
    columnFD = pd.MultiIndex.from_product(arraysFD, names=['region', 'category'])
    y_io = pd.DataFrame(np.random.random_sample((20, 8)), index, columnFD)

    dict_dataframes = {'F_io': F_io,  'C_io': C_io, 'STR_io': STR_io, 'y_io': y_io, 'A_io': A_io,
                       'io_database_name_and_version': io_database_name_and_version}
    return dict_dataframes

@pytest.fixture()
def random_parameters():
    listcountry = ['FR', 'CA', 'CH', 'ZA']
    listregions = ['GLO', 'RER']
    listcountry_per_regions = [listcountry, ['FR', 'CH']]

    replacements1 = {'CA-QC': 'CA'}

    reference_year_IO = 2011
    number_of_countries_IO = 4
    number_of_products_IO = 5

    dict_parameters = {'listcountry': listcountry, 'listregions': listregions,
                       'listcountry_per_regions': listcountry_per_regions, 'replacements1': replacements1,
                       'reference_year_IO': reference_year_IO, 'number_of_countries_IO': number_of_countries_IO,
                       'number_of_products_IO': number_of_products_IO}
    return dict_parameters
