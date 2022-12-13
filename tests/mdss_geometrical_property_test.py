import mdss.geometrical_property as g
import pytest


@pytest.mark.parametrize(
    "property_subclass , atom_selection",
    [
        (g.RMSDProperty, "name CA"),
        (
            g.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (g.RadiusOfGyrationProperty, "name CA"),
        (g.Angles, ["name CA", "name CA", "name CA"]),
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
    property_subclass, atom_selection, protein_data
):
    p_prop = property_subclass(protein_data, atom_selection)
    p_prop.calculate_property()
    assert len(p_prop.property_vector) == len(protein_data.frame_indices)


@pytest.mark.parametrize(
    "property_subclass , atom_selection",
    [
        (g.RMSDProperty, "name CA"),
        (
            g.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (g.RadiusOfGyrationProperty, "name CA"),
        (g.Angles, ["name CA", "name CA", "name CA"]),
    ],
)
def test_property_frame_indices_have_same_length(
    property_subclass, atom_selection, protein_data
):
    p_prop = property_subclass(protein_data, atom_selection)
    p_prop.calculate_property()
    assert len(p_prop.frame_indices) == len(protein_data.frame_indices)


@pytest.mark.parametrize(
    "property_subclass, atom_selection",
    [
        (g.DistanceBetweenAtoms, "name CA"),
        (g.Angles, "name CA"),
        (g.DihedralAngles, "name CA"),
    ],
)
def test_unexpected_atom_selection_raises_error(
    protein_data, property_subclass, atom_selection
):
    with pytest.raises(RuntimeError):
        property_subclass(protein_data, atom_selection)


@pytest.mark.parametrize(
    "property_subclass, atom_selection",
    [
        (g.RMSDProperty, "name CA"),
        (
            g.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (g.RadiusOfGyrationProperty, "name CA"),
        (g.Angles, ["name CA", "name CA", "name CA"]),
    ],
)
def test_property_statistics_method_have_been_called_once(
    property_subclass, atom_selection, protein_data, spy_method
):
    p_prop = property_subclass(protein_data, atom_selection)
    spy = spy_method(p_prop, "_property_statistics")
    p_prop.calculate_property()
    spy.assert_called_once()


@pytest.mark.parametrize(
    "property_subclass, atom_selection",
    [
        (g.RMSDProperty, "name CA"),
        (
            g.DistanceBetweenAtoms,
            ["name CA", "name CA"],
        ),
        (g.RadiusOfGyrationProperty, "name CA"),
        (g.Angles, ["name CA", "name CA", "name CA"]),
    ],
)
def test_discretize_vector_method_have_been_called_once(
    property_subclass, atom_selection, protein_data, spy_method
):
    p_prop = property_subclass(protein_data, atom_selection)
    spy = spy_method(p_prop, "discretize_vector")
    p_prop.calculate_property()
    spy.assert_called_once()
