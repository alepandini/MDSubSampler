import mdss_protein_data
import mdss_geometrical_property
from unittest.mock import MagicMock

"""       
Test 1: calculate_property() for RMSD
Lenth(self.property_vector) = length(self.protein_data.frame_indices)
length(frame_indices.append) = same
property_Statistics and discretized vector methods have been called once
pytest check how to test this
Test 2: calculate_property() for Distance between atoms

self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[0] called 2 times 
same as test 1 with length
Test 3: calculate_property() for Radius of Gyration

Test 4: calculate_property() for Angles

Test 5: calculate_property() for Dihedral angle phi

Test 6: calculate_property() for Dihedral angle psi


test that a function is called
assert that an error has been raised - give wrong arguments so we can have errors
"""


def test_RMSDProperty_vector_and_indices_have_same_lenght():

    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    p_prop = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    p_prop.calculate_property()
    assert len(p_prop.property_vector) == len(p_data.frame_indices)


def test_RMSDProperty_frame_indices_have_same_lenght():

    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    p_prop = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    p_prop.calculate_property()
    assert len(p_prop.frame_indices) == len(p_data.frame_indices)


# def test_RMSDProperty_property_statistics_method_have_been_called_once():

#     p_data = mdss_protein_data.ProteinData(
#         "data/user.xtc", "data/user.gro", config_parameters=None
#     )
#     p_prop = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
#     p_prop.calculate_property()
#     mock = mock.Mock()
#     mock.method()
#     mock.method.assert_called_once()


def test_DistanceBetweenAtoms_vector_and_indices_have_same_lenght():

    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    p_prop = mdss_geometrical_property.DistanceBetweenAtoms(
        p_data, atom_selection="name CA"
    )
    p_prop.calculate_property()
    assert len(p_prop.property_vector) == len(p_data.frame_indices)
