import random
import mdss_sampler
import mdss_dissimilarity
import mdss_protein_data
import mdss_geometrical_property
import mdss_property
import pytest
import os

# here = os.path.abspath(os.path.dirname(__file__))
# data_dir = os.path.join(here, "data")
# traj_file = os.path.join(data_dir, "user.xtc")
# top_file = os.path.join(data_dir, "user.gro")


def test_random_sampler_sample_has_expected_length():
    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    rmsd = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    rmsd.calculate_property()
    p_sampler = mdss_sampler.RandomSampler(rmsd, seed_number=1999)
    sampled_rmsd = p_sampler.sample(100)
    assert len(sampled_rmsd.frame_indices) == 100


def test_random_sampler_sample_returns_subset():
    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    rmsd = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    rmsd.calculate_property()
    frame_indices = p_data.frame_indices
    p_sampler = mdss_sampler.RandomSampler(rmsd, seed_number=1999)
    sampled_rmsd = p_sampler.sample(100)

    assert all(x in frame_indices for x in sampled_rmsd.frame_indices)


# def test_calculate_distance():
#     class FakeProperty:
#         def __init__(self):
#             self.avg_value = random.random()

#     property_1 = FakeProperty()
#     property_2 = FakeProperty()
#     distance = mdss_dissimilarity.Dissimilarity(property_1, property_2)

#     assert distance.distance == property_1.avg_value - property_2.avg_value


# @pytest.mark.parametrize(
#     "mocked_function_name, distance_subclass",
#     [
#         ("dictances.bhattacharyya", mdss_dissimilarity.BhattaCoefficient),
#         ("dictances.kullback_leibler", mdss_dissimilarity.KLDivergence),
#         ("dictances.pearson", mdss_dissimilarity.PearsonCoefficient),
#     ],
# )
# def test_bhattacharyya_distance(mocker, mocked_function_name, distance_subclass):
#     class FakeProperty:
#         def __init__(self):
#             self.property_vector_discretized = random.random()

#     mocked_function = mocker.patch(mocked_function_name)

#     property_1 = FakeProperty()
#     property_2 = FakeProperty()
#     distance_subclass(property_1, property_2)

#     mocked_function.assert_called_once_with(
#         property_1.property_vector_discretized, property_2.property_vector_discretized
#     )
