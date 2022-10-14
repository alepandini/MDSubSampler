import mdss_protein_data
import mdss_property

"""       
Test 2: discretize_vector()
check the values with an example of data, is it compliant with what I want to calculate
golden tests - maybe later

Test 5: write_property_vector()
Test if the property vector is written. Do we still have a property vector does this need to change?
give an out
Test 6: write_write_property_discretised_vector()
Same as test 5
"""
import tempfile
import os.path


def test_write_property_vector_method_has_written_the_file():
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfilepath = os.path.join(tmpdirname, "test.txt")
        p_data = mdss_protein_data.ProteinData(
            "data/user.xtc", "data/user.gro", config_parameters=None
        )
        p_prop = mdss_property.ProteinProperty(p_data, atom_selection="name CA")
        p_prop.write_property_vector(outfilepath)
        assert os.path.exists(outfilepath)
