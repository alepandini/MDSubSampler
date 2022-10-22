from audioop import rms
import mdss_dissimilarity
import mdss_sampler
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
    "dissimilarity_subclass, mocked_function",
    [
        (mdss_dissimilarity.Bhattacharya, "dictances.bhattacharyya"),
        (mdss_dissimilarity.Pearson, "dictances.pearson"),
    ],
)
def test_dissimilarity_measure(
    dissimilarity_subclass, mocked_function, mocker, rmsd_property
):
    sampler = mdss_sampler.RandomSampler(rmsd_property)
    sampled_property = sampler.sample(100)
    dissimilarity = dissimilarity_subclass(sampled_property, rmsd_property)
    function = mocker.patch(mocked_function)
    function.return_value = 0
    dissimilarity.calculate_dissimilarity()
    function.assert_called_once_with(
        dissimilarity.target_property.property_distribution_dict,
        dissimilarity.ref_property.property_distribution_dict,
    )
