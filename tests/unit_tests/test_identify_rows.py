import pytest
import pandas as pd
import numpy as np
import os
dir_name = os.path.dirname(os.path.realpath('pylcaio.py'))
import sys
sys.path.append(dir_name+'\\src')
import pylcaio

# implement quality checks with RoW left in io_geography, too many countries for IO resolution, GLO, etc.

def test_identify_rows(random_LCA_system, random_IO_system, random_parameters):
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
                                   listcountry_per_regions = random_parameters['listcountry_per_regions'],
                                   replacements1 = random_parameters['replacements1'],
                                   reference_year_IO = random_parameters['reference_year_IO'],
                                   number_of_countries_IO = random_parameters['number_of_countries_IO'],
                                   numer_of_products_IO = random_parameters['number_of_products_IO'])
    pylcaio_object.identify_rows()

    assert ['FR', 'CH', 'ZA'] in pylcaio_object.dictRoW.values()
    assert ['CA', 'ZA'] in pylcaio_object.dictRoW.values()
    assert len(pylcaio_object.dictRoW.values()) == 2
    assert type(pylcaio_object.dictRoW) == dict
    assert 'RoW' not in pylcaio_object.PRO_f.io_geography.tolist()

