import mdss_protein_data
import mdss_geometrical_property
import pytest


@pytest.mark.parametrize(
    "property_subclass , atom_selection",
    [
        (mdss_geometrical_property.RMSDProperty, "name CA"),
        (
            mdss_geometrical_property.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (mdss_geometrical_property.RadiusOfGyrationProperty, "name CA"),
        (mdss_geometrical_property.Angles, ["name CA", "name CA", "name CA"]),
        # (
        #     mdss_geometrical_property.DihedralAnglePhi,
        #     ["name CA", "name CA", "name CA", "name CA"],
        # ),
        # (
        #     mdss_geometrical_property.DihedralAnglePsi,
        #     ["name CA", "name CA", "name CA", "name CA"],
        # ),
    ],
)
def test_property_vector_and_indices_have_same_length(
    property_subclass, atom_selection, traj_file, top_file
):

    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    p_prop = property_subclass(p_data, atom_selection)
    p_prop.calculate_property()
    assert len(p_prop.property_vector) == len(p_data.frame_indices)


@pytest.mark.parametrize(
    "property_subclass , atom_selection",
    [
        (mdss_geometrical_property.RMSDProperty, "name CA"),
        (
            mdss_geometrical_property.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (mdss_geometrical_property.RadiusOfGyrationProperty, "name CA"),
        (mdss_geometrical_property.Angles, ["name CA", "name CA", "name CA"]),
    ],
)
def test_property_frame_indices_have_same_length(
    property_subclass, atom_selection, traj_file, top_file
):

    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    p_prop = property_subclass(p_data, atom_selection)
    p_prop.calculate_property()
    assert len(p_prop.frame_indices) == len(p_data.frame_indices)


@pytest.mark.parametrize(
    "property_subclass, atom_selection",
    [
        (mdss_geometrical_property.DistanceBetweenAtoms, "name CA"),
        (mdss_geometrical_property.Angles, "name CA"),
        (mdss_geometrical_property.DihedralAngles, "name CA"),
    ],
)
def test_unexpected_atom_selection_raises_error(
    traj_file, top_file, property_subclass, atom_selection
):
    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    with pytest.raises(RuntimeError):
        property_subclass(p_data, atom_selection)


@pytest.mark.parametrize(
    "property_subclass, atom_selection",
    [
        (mdss_geometrical_property.RMSDProperty, "name CA"),
        (
            mdss_geometrical_property.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (mdss_geometrical_property.RadiusOfGyrationProperty, "name CA"),
        (mdss_geometrical_property.Angles, ["name CA", "name CA", "name CA"]),
    ],
)
def test_property_statistics_method_have_been_called_once(
    property_subclass, atom_selection, traj_file, top_file, spy_method
):

    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    p_prop = property_subclass(p_data, atom_selection)
    spy = spy_method(p_prop, "_property_statistics")
    p_prop.calculate_property()
    spy.assert_called_once()


@pytest.mark.parametrize(
    "property_subclass, atom_selection",
    [
        (mdss_geometrical_property.RMSDProperty, "name CA"),
        (
            mdss_geometrical_property.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (mdss_geometrical_property.RadiusOfGyrationProperty, "name CA"),
        (mdss_geometrical_property.Angles, ["name CA", "name CA", "name CA"]),
    ],
)
def test_discretize_vector_method_have_been_called_once(
    property_subclass, atom_selection, traj_file, top_file, spy_method
):

    p_data = mdss_protein_data.ProteinData(traj_file, top_file, config_parameters=None)
    p_prop = property_subclass(p_data, atom_selection)
    spy = spy_method(p_prop, "discretize_vector")
    p_prop.calculate_property()
    spy.assert_called_once()
