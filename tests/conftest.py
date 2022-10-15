import pytest
import os
import os.path
import operator as op
import mdss_protein_data
import mdss_geometrical_property


@pytest.fixture
def tests_root():
    return os.path.abspath(os.path.dirname(__file__))


@pytest.fixture
def data_dir(tests_root):
    return os.path.join(tests_root, os.pardir, "data")


@pytest.fixture
def traj_file(data_dir):
    return os.path.join(data_dir, "user.xtc")


@pytest.fixture
def top_file(data_dir):
    return os.path.join(data_dir, "user.gro")


@pytest.fixture
def spy_method(mocker):
    """Python implementation for test spies.

    Works by mocking an object method and re-attaching the mocked
    method as a side-effect.
    """

    def _spy(object=None, method=None):
        return mocker.patch.object(
            object, method, side_effect=op.attrgetter(method)(object)
        )

    return _spy


@pytest.fixture
def protein_data(traj_file, top_file):
    return mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)


@pytest.fixture
def rmsd_property(protein_data):
    prop = mdss_geometrical_property.RMSDProperty(
        protein_data, atom_selection="name CA"
    )
    prop.calculate_property()
    return prop
