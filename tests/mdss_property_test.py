from os import rmdir
import tempfile
import os.path
import mdss.property as p
import mdss.sampler as s


def test_write_property_vector_method(rmsd_property):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfilepath = os.path.join(tmpdirname, "test.txt")
        rmsd_property.write_property_vector(outfilepath)
        assert os.path.exists(outfilepath)


def test_write_discretized_property_vector_method(rmsd_property):
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfilepath = os.path.join(tmpdirname, "test.txt")
        rmsd_property.write_discretized_property_vector(outfilepath)
        assert os.path.exists(outfilepath)


def test_get_samples_average(rmsd_property):
    p_sampler = s.ProteinSampler(rmsd_property)
    sampled_property = p.SampledProperty(
        rmsd_property,
        p_sampler.property_vector,
        p_sampler.frame_indices,
        p_sampler.samples_indices,
    )
    average_dict = sampled_property.get_samples_averages()
    assert len(average_dict) == len(set(p_sampler.samples_indices))


def test_get_samples_averages_method_called_in_get_averages(rmsd_property, spy_method):
    p_sampler = s.ProteinSampler(rmsd_property)
    sampled_property = p.SampledProperty(
        rmsd_property,
        p_sampler.property_vector,
        p_sampler.frame_indices,
        p_sampler.samples_indices,
    )
    spy = spy_method(sampled_property, "get_samples_averages")
    sampled_property.get_average()
    spy.assert_called_once()
