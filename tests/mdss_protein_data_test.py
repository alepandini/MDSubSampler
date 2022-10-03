import mdss_protein_data
import mdss_property
import mdss_geometrical_property
import mdss_parser as p
import mdss as m
import os.path
import numpy as np


here = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(here, "data")
traj_file = os.path.join(data_dir, "user.xtc")
top_file = os.path.join(data_dir, "user.gro")


def test_frames_of_trajectory_has_expected_length():

    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    frame_indices = p_data._frames_of_trajectory()
    assert len(frame_indices) == len(p_data.trajectory_data.trajectory)


def test_frame_selection_iterator_returns_selected_frames():
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    selection_of_frames = [0, 2, 4, slice(6, 10)]
    selected_frames = p_data.frame_selection_iterator(selection_of_frames)
    indices = [0, 2, 4, 6, 7, 8, 9]
    assert list(selected_frames) == [
        p_data.trajectory_data.trajectory[i] for i in indices
    ]


def test_frame_selection_indices_returns_selected_frames():
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    selection_of_frames = [0, 2, 4, slice(6, 10)]
    selected_frames = p_data.frame_selection_indices(selection_of_frames)
    assert list(selected_frames) == [0, 2, 4, 6, 7, 8, 9]


def test_add_property_is_linked_to_property_data():
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    p_prop = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    assert len(p_data.property_dict.keys()) == 1
    key = list(p_data.property_dict.keys())[0]
    assert key.startswith(p_prop.display_name)
    assert p_data.property_dict[key] is p_prop


"""       
Test 8: property_data_report(self) :
        Check that the values are inserted correctly in the dictionary
        (Make fake dictionary (key=rmsd, values=fake class that have min, max - with the testing framework)
        Enter fake input and check if expected output is there. 
Test 1: _frames_of_trajectory(self) : 
        Check that the expected output is returned (len(frames) = len(self.trajectory_data.trajectory))
        (Create a fake object and fake attributes then call function and see if we have the expected output)
"""
