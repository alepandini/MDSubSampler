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

with tempfile.TemporaryDirectory() as tmpdirname:
    outfilepath = os.path.join(tmpdirname, "test.txt")
    # create protein property
    # call write_property with outfilepath
    assert os.path.exists(outfilepath)
