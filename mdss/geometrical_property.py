from mdss.property import ProteinProperty
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import dihedrals
from mdss.log_setup import log
import numpy as np


class RMSD(ProteinProperty):
    """
    Represents RMSD property class
    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    atom_selection: str
        Atom selection for property calculation
    fit : performs RMSD superposition if True
    """

    display_name = "RMSD"

    def __init__(self, protein_data, atom_selection="name CA", fit=False):
        self.fit = fit
        super().__init__(protein_data, atom_selection)

    def calculate_property(self, frame_index=None):
        """
        Calculates RMSD property of frame list for a selection of atoms
        """

        if self.set_reference_coordinates(frame_index):
            for frame in self.protein_data.frame_indices:
                """
                Go through the trajectory and for each frame I compare with my reference frame
                """
                self.protein_data.trajectory_data.trajectory[frame]

                self.property_vector.append(
                    rms.rmsd(
                        self.protein_data.trajectory_data.select_atoms(
                            self.atom_selection
                        ).positions,
                        self.ref_coordinates,
                        center=self.fit,
                        superposition=self.fit,
                    )
                )
                self.frame_indices.append(frame)

            self._property_statistics()
            self.discretize_vector()
        else:
            print("Property cannot be calculated without associated protein data")
            log.warning(
                "{:12s} Property cannot be calculated without associated protein data".format(
                    "STEPS"
                )
            )


class DistanceBetweenAtoms(ProteinProperty):
    """
    Represents Distance between atoms property class

    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    atom_selection: list
        List with 2 atoms or 2 group of atoms
    """

    display_name = "Distance_between_atoms"

    def __init__(self, protein_data, atom_selection):
        if not isinstance(atom_selection, list) or len(atom_selection) != 2:
            log.error(
                "{:18s} Expecting atom_selection to be a list of 2 selections in DistanceBetweenAtoms class".format(
                    "INPUT"
                )
            )
            raise RuntimeError("Expecting atom_selection to be a list of 2 selections")

        super().__init__(protein_data, atom_selection)

    def calculate_property(self, frame_index=None):
        """
        Calculates distance between two given set of atoms
        """
        if self.set_reference_coordinates(frame_index):
            atom_selection_1 = self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[0]
            )
            atom_selection_2 = self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[1]
            )

            for frame in self.protein_data.frame_indices:
                """
                Go through the trajectory and for each frame the distance between the given
                atoms is calculated
                """
                self.protein_data.trajectory_data.trajectory[frame]
                dist = distances.distance_array(
                    atom_selection_1.positions[0], atom_selection_2.positions[1]
                )
                self.property_vector.append(dist)
                self.frame_indices.append(frame)

            self._property_statistics()
            self.discretize_vector()
        else:
            print("Property cannot be calculated without associated protein data")
            log.warning(
                "{:12s} Property cannot be calculated without associated protein data".format(
                    "STEPS"
                )
            )


class RadiusOfGyrationProperty(ProteinProperty):
    """
    Represents Radius of Gyration property class
    """

    display_name = "Radius of Gyration"

    def calculate_property(self):
        """
        Calculates radius of gyration of atoms in each frame
        """
        self.time = []
        for frame in self.protein_data.frame_indices:
            """
            Go through the trajectory and for the atoms of each frame the rog is calculated
            """
            self.protein_data.trajectory_data.trajectory[frame]
            self.time.append(self.protein_data.trajectory_data.trajectory.time)
            self.property_vector.append(
                self.protein_data.trajectory_data.select_atoms(
                    self.atom_selection
                ).radius_of_gyration()
            )
            self.frame_indices.append(frame)

        self._property_statistics()
        self.discretize_vector()


class Angles(ProteinProperty):
    """
    Represents Angles property class

    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    atom_selection: list
       List with selection of 3 atoms
    """

    display_name = "Angle between 3 atoms"

    def __init__(self, protein_data, atom_selection):
        if not isinstance(atom_selection, list) or len(atom_selection) != 3:
            log.error(
                "{:18s} Expecting atom_selection to be a list of 3 selections in Angles class".format(
                    "INPUT"
                )
            )
            raise RuntimeError("Expecting atom_selection to be a list of 3 selections")

        super().__init__(protein_data, atom_selection)

    def calculate_property(self, frame_index=None):
        """
        Calculates angle between three given atoms
        """
        if self.set_reference_coordinates(frame_index):
            atom_selection_1 = self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[0]
            )
            atom_selection_2 = self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[1]
            )
            atom_selection_3 = self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[2]
            )

            for frame in self.protein_data.frame_indices:
                """
                Go through the trajectory and for each frame the angle between the given
                atoms is calculated
                """

                self.protein_data.trajectory_data.trajectory[frame]
                atom_1_coordinates = atom_selection_1.positions[0]
                atom_2_coordinates = atom_selection_2.positions[1]
                atom_3_coordinates = atom_selection_3.positions[2]

                atoms_2_1 = atom_1_coordinates - atom_2_coordinates
                atoms_2_3 = atom_3_coordinates - atom_2_coordinates

                cosine_angle = np.dot(atoms_2_1, atoms_2_3) / (
                    np.linalg.norm(atoms_2_1) * np.linalg.norm(atoms_2_3)
                )
                angle = np.arccos(cosine_angle)
                angle_degrees = np.degrees(angle)
                self.property_vector.append(angle_degrees)
                self.frame_indices.append(frame)

            self._property_statistics()
            self.discretize_vector()
        else:
            print("Property cannot be calculated without associated protein data")
            log.warning(
                "{:12s} Property cannot be calculated without associated protein data".format(
                    "STEPS"
                )
            )


class DihedralAngles(ProteinProperty):
    """
    Represents dihedral angles property class

    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    atom_selection: list
        List with selection of 4 atoms
    """

    def __init__(self, protein_data, atom_selection):
        if not isinstance(atom_selection, list) or len(atom_selection) != 4:
            log.error(
                "{:18s} Expecting atom_selection to be a list of 4 selections in DihedralAngles class".format(
                    "INPUT"
                )
            )
            raise RuntimeError("Expecting atom_selection to be a list of 4 selections")

        super().__init__(protein_data, atom_selection)


class DihedralAnglePhi(DihedralAngles):
    """
    Represents Dihedral angle phi property class
    """

    display_name = "Dihedral Angle phi between 4 selected atoms"

    def calculate_property(self, frame_index=None):
        if self.set_reference_coordinates(frame_index):
            u = self.protein_data.trajectory_data
            for frame in self.protein_data.frame_indices:
                phi_ags = [res.phi_selection() for res in u.residues]
                phi_ags = [phi for phi in phi_ags if phi is not None]
                dihs = dihedrals.Dihedral(phi_ags).run()
                self.property_vector.append(dihs.results.angles)
                self.frame_indices.append(frame)

            # self._property_statistics()
            # self.discretize_vector()
        else:
            print("Property cannot be calculated without associated protein data")
            log.warning(
                "{:12s} Property cannot be calculated without associated protein data".format(
                    "STEPS"
                )
            )


class DihedralAnglePsi(DihedralAngles):
    """
    Represents Dihedral angle psi property class
    """

    display_name = "Dihedral Angle psi between 4 selected atoms"

    def calculate_property(self, frame_index=None):
        if self.set_reference_coordinates(frame_index):
            u = self.protein_data.trajectory_data
            for frame in self.protein_data.frame_indices:
                psi_ags = [res.psi_selection() for res in u.residues]
                psi_ags = [psi for psi in psi_ags if psi is not None]
                dihs = dihedrals.Dihedral(psi_ags).run()
                self.property_vector.append(dihs.results.angles)
                self.frame_indices.append(frame)

            # self._property_statistics()
            # self.discretize_vector()
        else:
            print("Property cannot be calculated without associated protein data")
            log.warning(
                "{:12s} Property cannot be calculated without associated protein data".format(
                    "STEPS"
                )
            )
