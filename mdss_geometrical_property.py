from mdss_property import ProteinProperty
from mdss_property import SampledProperty
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import dihedrals
from mdss_logging import logging as log
import numpy as np


class RMSDProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the RMSD value for each frame in the
    protein trajectory
    """

    display_name = "RMSD"

    def calculate_property(self):
        """
        Method that calculates the RMDS proporty of a given frame list for a given
        selection of atoms
        """
        if self.set_reference_coordinates():
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
                    )
                )
                self.frame_indices.append(frame)

            self._property_statistics()
            self.discretize_vector()
        else:
            print("Property cannot be calculated without associated protein data")
            log.warning(
                "{:18s} Property cannot be calculated without associated protein data".format(
                    "WARNING"
                )
            )


class DistanceBetweenAtoms(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the distance value between two
    atoms or two group of atoms for each frame in the protein trajectory

    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class

    atom_selection: list
        A list with selection of atoms for distance calculation between them
    """

    display_name = "Distance_between_atoms"

    def __init__(self, protein_data, atom_selection):
        if not isinstance(atom_selection, list) or len(atom_selection) != 2:
            log.error(
                "{:18s} Expecting atom_selection to be a list of 2 selections in DistanceBetweenAtoms class".format(
                    "ERROR"
                )
            )
            raise RuntimeError("Expecting atom_selection to be a list of 2 selections")

        super().__init__(protein_data, atom_selection)

    def calculate_property(self):
        """
        Method that calculates the distance between two given set of atoms
        """
        if self.set_reference_coordinates():
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
                "{:18s} Property cannot be calculated without associated protein data".format(
                    "WARNING"
                )
            )


class RadiusOfGyrationProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the Radius of Gyration value
    for each frame in the protein trajectory
    """

    display_name = "Radius of Gyration"

    def calculate_property(self):
        """
        Method that calculates the radius of gyration of the atoms in each frame
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
    A Subclass of ProteinProperty class that calculates the angles between 3 selected atoms
    in the protein structure

    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class

    atom_selection: list
        A list with selection of atoms for calculation of angle they form
    """

    display_name = "Angle between 3 atoms"

    def __init__(self, protein_data, atom_selection):
        if not isinstance(atom_selection, list) or len(atom_selection) != 3:
            log.error(
                "{:18s} Expecting atom_selection to be a list of 3 selections in Angles class".format(
                    "ERROR"
                )
            )
            raise RuntimeError("Expecting atom_selection to be a list of 3 selections")

        super().__init__(protein_data, atom_selection)

    def calculate_property(self):
        """
        Method that calculates the angle between three given atoms
        """
        if self.set_reference_coordinates():
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
                "{:18s} Property cannot be calculated without associated protein data".format(
                    "WARNING"
                )
            )


class DihedralAngles(ProteinProperty):
    """
    A Subclass of ProteinProperty class that calculates the dihedral angle
    between 4 selected atoms in the protein structure

    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class

    atom_selection: list
        A list with selection of atoms for calculation of dihedral angle they form
    """

    def __init__(self, protein_data, atom_selection):
        if not isinstance(atom_selection, list) or len(atom_selection) != 4:
            log.error(
                "{:18s} Expecting atom_selection to be a list of 4 selections in DihedralAngles class".format(
                    "ERROR"
                )
            )
            raise RuntimeError("Expecting atom_selection to be a list of 4 selections")

        super().__init__(protein_data, atom_selection)


class DihedralAnglePhi(DihedralAngles):
    """
    A Subclass of DihedralAngles class that calculates the dihedral angle phi
    between 4 selected atoms in the protein structure
    """

    display_name = "Dihedral Angle phi between 4 selected atoms"

    def calculate_property(self):
        if self.set_reference_coordinates():
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
                "{:18s} Property cannot be calculated without associated protein data".format(
                    "WARNING"
                )
            )


class DihedralAnglePsi(DihedralAngles):
    """
    A Subclass of DihedralAngles class that calculates the dihedral angle psi
    between 4 selected atoms in the protein structure
    """

    display_name = "Dihedral Angle psi between 4 selected atoms"

    def calculate_property(self):
        if self.set_reference_coordinates():
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
                "{:18s} Property cannot be calculated without associated protein data".format(
                    "WARNING"
                )
            )
