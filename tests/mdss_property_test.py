import tempfile
import os.path

"""       
Test discretize_vector()
golden tests - maybe later

Tests methods of SampledProperty later

"""


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
