import pytest
import pandas as pd
import numpy as np
import os
dir_name = os.path.dirname(os.path.realpath('pylcaio.py'))
import sys
sys.path.append(dir_name+'\\src')
import pylcaio


def test_hybridize(random_LCA_system, random_IO_system, random_parameters):
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
                                   listmarket = random_LCA_system['listmarket'],
                                   number_of_RoW_IO = random_parameters['number_of_RoW_IO'],
                                   io_categories=random_parameters['io_categories'],
                                   STAM_table=random_parameters['STAM_table'],
                                   categories_same_functionality=random_parameters['categories_same_functionality'])
    pylcaio_object.hybridize('STAM')

    assert pylcaio_object.A_io_f.sum().sum() != 0
    assert not pylcaio_object.A_io_f.empty
    assert (pylcaio_object.A_io_f.index == pylcaio_object.A_io.index).all()
    assert (pylcaio_object.A_io_f.columns == pylcaio_object.A_ff.columns).all()

    # all the sector2, 3 and 5 rows should be zero
    assert pylcaio_object.A_io_f.loc[[i for i in pylcaio_object.A_io_f.index if i[1] == 'sector2']].sum().sum() == 0
    assert pylcaio_object.A_io_f.loc[[i for i in pylcaio_object.A_io_f.index if i[1] == 'sector3']].sum().sum() == 0
    assert pylcaio_object.A_io_f.loc[[i for i in pylcaio_object.A_io_f.index if i[1] == 'sector5']].sum().sum() == 0

    # the sector4 row should not have ANY zero value
    assert pylcaio_object.A_io_f.loc[[i for i in pylcaio_object.A_io_f.index if i[1] == 'sector4']][
        pylcaio_object.A_io_f.loc[[i for i in pylcaio_object.A_io_f.index if i[1] == 'sector4']]==0].isnull().all().all()

    # no zero in row sector1 columns 0, 2 and 4
    assert pylcaio_object.A_io_f.loc[[i for i in pylcaio_object.A_io_f.index if i[1] == 'sector1'],
                                     [pylcaio_object.A_io_f.columns[0],
                                      pylcaio_object.A_io_f.columns[2],
                                      pylcaio_object.A_io_f.columns[4]]][
        pylcaio_object.A_io_f.loc[
            [i for i in pylcaio_object.A_io_f.index if i[1] == 'sector1'],[pylcaio_object.A_io_f.columns[0],
                                                                           pylcaio_object.A_io_f.columns[2],
                                                                           pylcaio_object.A_io_f.columns[4]]]==0].isnull().all().all()
    # only zeros in row sector1 columns 1 and 3
    assert pylcaio_object.A_io_f.loc[[i for i in pylcaio_object.A_io_f.index if i[1] == 'sector1'],[
        pylcaio_object.A_io_f.columns[1],pylcaio_object.A_io_f.columns[3]]].sum().sum() == 0

    # check that if uncorrected is different from corrected, then it is equal to zero (with STAM)
    assert pylcaio_object.A_io_f.loc[
        [i for i in pylcaio_object.A_io_f_uncorrected.index
         if pylcaio_object.A_io_f_uncorrected.loc[
             i,pylcaio_object.A_io_f_uncorrected.columns[0]] != pylcaio_object.A_io_f.loc[
             i,pylcaio_object.A_io_f_uncorrected.columns[0]]]].sum().sum() == 0