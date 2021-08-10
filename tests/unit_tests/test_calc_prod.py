import pytest
import pandas as pd
import numpy as np
import os
dir_name = os.path.dirname(os.path.realpath('pylcaio.py'))
import sys
sys.path.append(dir_name+'\\src')
import pylcaio


def test_calc_proc(random_LCA_system, random_IO_system, random_parameters):
    pylcaio_object = pylcaio.LCAIO(PRO_f = random_LCA_system['PRO_f'], A_ff = random_LCA_system['A_ff'],
                                   F_f = random_LCA_system['F_f'], C_f = random_LCA_system['C_f'],
                                   STR_f = random_LCA_system['STR_f'], y_f = random_LCA_system['y_f'],
                                   lca_database_name_and_version = 'ecoinvent3.5',
                                   A_io = random_IO_system['A_io'], F_io = random_IO_system['F_io'],
                                   C_io = random_IO_system['C_io'], STR_io = random_IO_system['STR_io'],
                                   y_io = random_IO_system['y_io'],
                                   io_database_name_and_version = 'exiobase3',
                                   listcountry = random_parameters['listcountry'],
                                   listregions = random_parameters['listregions'],
                                   countries_per_regions = random_parameters['countries_per_regions'],
                                   replacements1 = random_parameters['replacements1'],
                                   reference_year_IO = random_parameters['reference_year_IO'],
                                   number_of_countries_IO = random_parameters['number_of_countries_IO'],
                                   number_of_products_IO = random_parameters['number_of_products_IO'],
                                   list_to_hyb = random_LCA_system['list_to_hyb'],
                                   list_not_hyb = random_LCA_system['list_not_to_hyb'],
                                   listmarket = random_LCA_system['listmarket'])
    pylcaio_object.identify_rows()
    pylcaio_object.calc_productions()

    # check prod not empty
    assert len(pylcaio_object.total_prod_country) != 0
    assert len(pylcaio_object.total_prod_region) != 0
    assert len(pylcaio_object.total_prod_RoW) != 0
    # check no string in prod
    assert type(pylcaio_object.total_prod_country.sum().sum()) == np.float64
    assert type(pylcaio_object.total_prod_region.sum().sum()) == np.float64
    assert type(pylcaio_object.total_prod_RoW.sum().sum()) == np.float64
    # check non NaN in prod
    assert not pylcaio_object.total_prod_country.isnull().any()[0]
    assert not pylcaio_object.total_prod_region.isnull().any()[0]
    assert not pylcaio_object.total_prod_RoW.isnull().any()[0]
    # check if prods are in a dataframe format
    assert type(pylcaio_object.total_prod_country) == pd.core.frame.DataFrame
    assert type(pylcaio_object.total_prod_region) == pd.core.frame.DataFrame
    assert type(pylcaio_object.total_prod_RoW) == pd.core.frame.DataFrame
    # check no prod is negative
    assert pylcaio_object.total_prod_country[pylcaio_object.total_prod_country<0].sum()[0] == 0
    assert pylcaio_object.total_prod_region[pylcaio_object.total_prod_region<0].sum().sum() == 0
    assert pylcaio_object.total_prod_RoW[pylcaio_object.total_prod_RoW<0].sum().sum() == 0
    # check that indexes are the same
    assert (pylcaio_object.A_io.index == pylcaio_object.total_prod_country.index).all()
    assert (pylcaio_object.A_io.index.levels[1] == pylcaio_object.total_prod_region.index).all()
    assert (pylcaio_object.A_io.index.levels[1] == pylcaio_object.total_prod_RoW.index).all()
