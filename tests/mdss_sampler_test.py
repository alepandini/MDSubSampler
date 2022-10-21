import mdss_sampler
import random
import math
import pytest


@pytest.mark.parametrize(
    "sampler_subclass, kwargs, expected",
    [
        (mdss_sampler.RandomSampler, {}, 100),
        (mdss_sampler.UniformSampler, {"strata_number": 10}, 99),
        (mdss_sampler.BootstrappingSampler, {"number_of_iterations": 10}, 1000),
    ],
)
def test_sampler_sample_has_expected_length(
    rmsd_property, sampler_subclass, kwargs, expected
):
    p_sampler = sampler_subclass(rmsd_property, **kwargs)
    sampled_prop = p_sampler.sample(100)
    assert len(sampled_prop.frame_indices) == expected


def test_stratified_sampler_sample_has_expected_length(rmsd_property):
    s_vector = [
        round((v + random.random()) * 13) for v in rmsd_property.property_vector
    ]
    p_sampler = mdss_sampler.StratifiedSampler(rmsd_property, strata_vector=s_vector)
    sampled_prop = p_sampler.sample(100)
    assert math.isclose(len(sampled_prop.frame_indices), 100, rel_tol=3)


@pytest.mark.parametrize(
    "sampler_subclass, kwargs, expected",
    [
        (mdss_sampler.RandomSampler, {}, 100),
        (mdss_sampler.UniformSampler, {"strata_number": 10}, 99),
        (mdss_sampler.BootstrappingSampler, {"number_of_iterations": 10}, 1000),
    ],
)
def test_sampler_sample_returns_subset(
    protein_data, rmsd_property, sampler_subclass, kwargs, expected
):
    p_sampler = sampler_subclass(rmsd_property, **kwargs)
    sampled_prop = p_sampler.sample(100)
    frame_indices = protein_data.frame_indices
    assert all(x in frame_indices for x in sampled_prop.frame_indices)


def test_stratified_sampler_sample_returns_subset(protein_data, rmsd_property):
    frame_indices = protein_data.frame_indices
    s_vector = [
        round((v + random.random()) * 13) for v in rmsd_property.property_vector
    ]
    p_sampler = mdss_sampler.StratifiedSampler(rmsd_property, s_vector)
    sampled_prop = p_sampler.sample(100)
    assert all(x in frame_indices for x in sampled_prop.frame_indices)


@pytest.mark.parametrize("size", [100, "100", "30%", "30.8%"])
def test_sample_calls_implementation(rmsd_property, spy_method, size):
    p_sampler = mdss_sampler.ProteinSampler(rmsd_property)
    spy = spy_method(p_sampler, "_sample")
    p_sampler.sample(size)
    spy.assert_called_once_with(
        mdss_sampler.convert_size(size, rmsd_property.protein_data.n_frames)
    )


@pytest.mark.parametrize(
    "given, expected, n_frames",
    [
        (100, 100, 1),
        ("100", 100, 1),
        ("30%", 30, 100),
        ("30%", 30, 101),
        ("30.8%", 31, 100),
    ],
)
def test_convert_size(given, expected, n_frames):
    assert mdss_sampler.convert_size(given, n_frames) == expected
