import pytest


# check that these potential names of database pass
@pytest.mark.parametrize('test_names_to_succeed', ['ecoinvent3.5', 'ecoinvent3.3', 'ecoinvent3.4', 'Eora45',
                                                   'exiobase3', 'exiobase2', '0.8'])
def test_extract_version_from_name_success(test_names_to_succeed):

    version = extract_version_from_name(test_names_to_succeed)
    try:
        assert int(version)
    except ValueError:
        assert float(version)


# function to test (from pylcaio)
def extract_version_from_name(name_database):

    for i in range(0, len(name_database)):
        try:
            if type(int(name_database[i])) == int:
                return name_database[i:]
        except ValueError:
            pass
