def test_frames_of_trajectory_has_expected_length(protein_data):
    frame_indices = protein_data._frames_of_trajectory()
    assert len(frame_indices) == len(protein_data.trajectory_data.trajectory)


def test_frame_selection_iterator_returns_selected_frames(protein_data):
    selection_of_frames = [0, 2, 4, slice(6, 10)]
    selected_frames = protein_data.frame_selection_iterator(selection_of_frames)
    indices = [0, 2, 4, 6, 7, 8, 9]
    assert list(selected_frames) == [
        protein_data.trajectory_data.trajectory[i] for i in indices
    ]


def test_frame_selection_indices_returns_selected_frames(protein_data):
    selection_of_frames = [0, 2, 4, slice(6, 10)]
    selected_frames = protein_data.frame_selection_indices(selection_of_frames)
    assert list(selected_frames) == [0, 2, 4, 6, 7, 8, 9]


def test_add_property_is_linked_to_property_data(protein_data, rmsd_property):
    assert len(protein_data.property_dict.keys()) == 1
    key = list(protein_data.property_dict.keys())[0]
    assert key.startswith(rmsd_property.display_name)
    assert protein_data.property_dict[key] is rmsd_property


def test_property_data_report(protein_data, rmsd_property):
    data_report = protein_data.property_data_report()
    assert len(data_report.keys()) == 1
    key = list(data_report.keys())[0]
    assert key.startswith(rmsd_property.display_name)
    assert data_report[key] == {
        "min": rmsd_property.min_value,
        "max": rmsd_property.max_value,
        "atom_selection": rmsd_property.atom_selection,
        "name": rmsd_property.display_name,
    }
