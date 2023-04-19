import mdss.sampler as s
import random
import math
import pytest
import tempfile
import os


@pytest.mark.parametrize(
    "sampler_subclass, kwargs, expected",
    [
        (s.RandomSampler, {}, 100),
        (s.UniformSampler, {"strata_number": 10}, 100),
        (s.WeightedSampler, {"weights_vector": None}, 100),
        (s.BootstrappingSampler, {"number_of_iterations": 10}, 1000),
    ],
)
def test_sampler_sample_has_expected_length(
    rmsd_property, sampler_subclass, kwargs, expected
):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfilepath = os.path.join(tmpdirname, "test.txt")
        file_prefix = "test"
        p_sampler = sampler_subclass(
            rmsd_property,
            rmsd_property.protein_data,
            outfilepath,
            file_prefix,
            **kwargs,
        )
        sampled_prop = p_sampler.sample(100)
        assert len(sampled_prop.frame_indices) == expected


def test_stratified_sampler_sample_has_expected_length(rmsd_property):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfilepath = os.path.join(tmpdirname, "test.txt")
        file_prefix = "test"
        s_vector = [
            round((v + random.random()) * 13) for v in rmsd_property.property_vector
        ]
        p_sampler = s.StratifiedSampler(
            rmsd_property,
            rmsd_property.protein_data,
            outfilepath,
            file_prefix,
            strata_vector=s_vector,
        )
        sampled_prop = p_sampler.sample(100)
        assert math.isclose(len(sampled_prop.frame_indices), 100, rel_tol=3)


@pytest.mark.parametrize(
    "sampler_subclass, kwargs, expected",
    [
        (s.RandomSampler, {}, 100),
        (s.UniformSampler, {"strata_number": 10}, 100),
        (s.WeightedSampler, {"weights_vector": None}, 100),
        (s.BootstrappingSampler, {"number_of_iterations": 10}, 1000),
    ],
)
def test_sampler_sample_returns_subset(
    protein_data, rmsd_property, sampler_subclass, kwargs, expected
):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfilepath = os.path.join(tmpdirname, "test.txt")
        file_prefix = "test"
        p_sampler = sampler_subclass(
            rmsd_property,
            rmsd_property.protein_data,
            outfilepath,
            file_prefix,
            **kwargs,
        )
        sampled_prop = p_sampler.sample(100)
        frame_indices = protein_data.frame_indices
        assert all(x in frame_indices for x in sampled_prop.frame_indices)


def test_stratified_sampler_sample_returns_subset(protein_data, rmsd_property):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfilepath = os.path.join(tmpdirname, "test.txt")
        file_prefix = "test"
        frame_indices = protein_data.frame_indices
        s_vector = [
            round((v + random.random()) * 13) for v in rmsd_property.property_vector
        ]
        p_sampler = s.StratifiedSampler(
            rmsd_property,
            rmsd_property.protein_data,
            outfilepath,
            file_prefix,
            strata_vector=s_vector,
        )
        sampled_prop = p_sampler.sample(100)
        assert all(x in frame_indices for x in sampled_prop.frame_indices)
