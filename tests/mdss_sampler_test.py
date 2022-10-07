import mdss_sampler
import mdss_protein_data
import mdss_geometrical_property
import random
import math

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


def test_stratified_sampler_sample_has_expected_length():
    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    rmsd = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    rmsd.calculate_property()
    s_vector = [round((v + random.random()) * 13) for v in rmsd.property_vector]
    p_sampler = mdss_sampler.StratifiedSampler(rmsd, s_vector)
    sampled_rmsd = p_sampler.sample(100)
    assert math.isclose(len(sampled_rmsd.frame_indices), 100, rel_tol=3)


def test_stratified_sampler_sample_returns_subset():
    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    rmsd = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    rmsd.calculate_property()
    frame_indices = p_data.frame_indices
    s_vector = [round((v + random.random()) * 13) for v in rmsd.property_vector]
    p_sampler = mdss_sampler.StratifiedSampler(rmsd, s_vector)
    sampled_rmsd = p_sampler.sample(100)
    assert all(x in frame_indices for x in sampled_rmsd.frame_indices)


def test_uniform_sampler_sample_has_expected_length():
    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    rmsd = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    rmsd.calculate_property()
    strata = 50
    p_sampler = mdss_sampler.UniformSampler(rmsd, strata)
    sampled_rmsd = p_sampler.sample(100)
    assert math.isclose(len(sampled_rmsd.frame_indices), 100, rel_tol=3)


def test_uniform_sampler_sample_returns_subset():
    p_data = mdss_protein_data.ProteinData(
        "data/user.xtc", "data/user.gro", config_parameters=None
    )
    rmsd = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    rmsd.calculate_property()
    strata = 50
    frame_indices = p_data.frame_indices
    p_sampler = mdss_sampler.UniformSampler(rmsd, strata)
    sampled_rmsd = p_sampler.sample(100)
    assert all(x in frame_indices for x in sampled_rmsd.frame_indices)
