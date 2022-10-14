import mdss_protein_data
import mdss_geometrical_property


def test_frames_of_trajectory_has_expected_length(traj_file, top_file):

    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    frame_indices = p_data._frames_of_trajectory()
    assert len(frame_indices) == len(p_data.trajectory_data.trajectory)


def test_frame_selection_iterator_returns_selected_frames(traj_file, top_file):
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    selection_of_frames = [0, 2, 4, slice(6, 10)]
    selected_frames = p_data.frame_selection_iterator(selection_of_frames)
    indices = [0, 2, 4, 6, 7, 8, 9]
    assert list(selected_frames) == [
        p_data.trajectory_data.trajectory[i] for i in indices
    ]


def test_frame_selection_indices_returns_selected_frames(traj_file, top_file):
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    selection_of_frames = [0, 2, 4, slice(6, 10)]
    selected_frames = p_data.frame_selection_indices(selection_of_frames)
    assert list(selected_frames) == [0, 2, 4, 6, 7, 8, 9]


def test_add_property_is_linked_to_property_data(traj_file, top_file):
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    p_prop = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    assert len(p_data.property_dict.keys()) == 1
    key = list(p_data.property_dict.keys())[0]
    assert key.startswith(p_prop.display_name)
    assert p_data.property_dict[key] is p_prop


def test_property_data_report(traj_file, top_file):
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    p_prop = mdss_geometrical_property.RMSDProperty(p_data, atom_selection="name CA")
    p_prop.calculate_property()
    data_report = p_data.property_data_report()
    assert len(data_report.keys()) == 1
    key = list(data_report.keys())[0]
    assert key.startswith(p_prop.display_name)
    assert data_report[key] == {
        "min": p_prop.min_value,
        "max": p_prop.max_value,
        "atom_selection": p_prop.atom_selection,
        "name": p_prop.display_name,
    }
