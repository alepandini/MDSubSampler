from audioop import rms
import mdss_dissimilarity
import mdss_sampler
import random
import pytest


def test_calculate_dissimilarity(rmsd_property):
    sampler = mdss_sampler.RandomSampler(rmsd_property)
    sampled_property = sampler.sample(100)
    dissimilarity = mdss_dissimilarity.Dissimilarity(sampled_property, rmsd_property)
    assert (
        dissimilarity.calculate_dissimilarity()
        == sampled_property.avg_value - rmsd_property.avg_value
    )


@pytest.mark.parametrize(
    "dissimilarity_subclass",
    [
        (mdss_dissimilarity.Bhattacharya),
        (mdss_dissimilarity.KullbackLeibler),
        (mdss_dissimilarity.Pearson),
    ],
)
def test_dissimilarity_measure(dissimilarity_subclass, rmsd_property, spy_method):
    sampler = mdss_sampler.RandomSampler(rmsd_property)
    sampled_property = sampler.sample(100)
    dissimilarity = dissimilarity_subclass(sampled_property, rmsd_property)
    spy = spy_method(dissimilarity, "calculate_dissimilarity")
    dissimilarity.calculate_dissimilarity()
    # spy.assert_called_once_with(
    #     sampled_property.property_vector,
    #     rmsd_property.property_vector,
    # )
    spy.assert_called_once()


# def test_calculate_dissimilarity():
#     class FakeProperty:
#         def __init__(self):
#             self.avg_value = random.random()

#     target_property = FakeProperty()
#     ref_property = FakeProperty()
#     dissimilarity = mdss_dissimilarity.Dissimilarity(target_property, ref_property)

#     assert (
#         dissimilarity.calculate_dissimilarity
#         == target_property.avg_value - ref_property.avg_value
#     )


# @pytest.mark.parametrize(
#     "mocked_function_name, distance_subclass",
#     [
#         ("dictances.bhattacharyya", mdss_dissimilarity.Bhattacharya),
#         ("dictances.kullback_leibler", mdss_dissimilarity.KullbackLeibler),
#         ("dictances.pearson", mdss_dissimilarity.Pearson),
#     ],
# )
# def test_bhattacharyya_distance(mocker, mocked_function_name, distance_subclass):
#     class FakeProperty:
#         def __init__(self):
#             self.property_distribution_dict = random.random()

#     mocked_function = mocker.patch(mocked_function_name)

#     target_property = FakeProperty()
#     ref_property = FakeProperty()
#     distance_subclass(target_property, ref_property)

#     mocked_function.assert_called_once_with(
#         target_property.property_distribution_dict,
#         ref_property.property_distribution_dict,
#     )
