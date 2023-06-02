"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Namir Oues
    
    This file is part of the MDSubSampler software 
    (https://github.com/alepandini/MDSubSampler).
    Copyright (c) 2023 Namir Oues and Alessandro Pandini.

    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
from mdss.property import ProteinProperty
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import dihedrals
from mdss.log_setup import log
import numpy as np


class RMSD(ProteinProperty):
    """
    Subclass of ProteinProperty class representing RMSD geometric property.

    Parameters
    ----------
    protein_data   : ProteinData
                     An instance of the ProteinData class representing the protein data.
    atom_selection : str
                     Atom selection for property calculation.
    fit            : boolean
                     If True, performs superposition of all structures to the
                     referene structure before RMSD calculation.
    """

    display_name = "RMSD"

    def __init__(self, protein_data, atom_selection="name CA", fit=False):
        self.fit = fit
        super().__init__(protein_data, atom_selection)

    def calculate_property(self, frame_index=None):
        """
        Calculate RMSD property for all trajectory frames for a selection of atoms.

        Parameters
        ----------
        frame_index : int
                      Reference structure (i.e. frame) from inputed protein trajectory.
        """

        if self.set_reference_coordinates(frame_index):
            for frame in self.protein_data.frame_indices:
                # Goes through the trajectory and for each frame I compare with my reference frame
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
    Subclass of ProteinProperty class representing distance geometric property between
    2 atoms or 2 group of atoms.

    Parameters
    ----------
    protein_data   : ProteinData
                     An instance of the ProteinData class representing the protein data.
    atom_selection : list
                     List of atom selection with 2 atoms or 2 group of atoms.
    """

    display_name = "DBA"

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
        Calculate distance between 2 given set of atoms for all trajectory frames.

        Parameters
        ----------
        frame_index : int
                      Reference structure (i.e. frame) from inputed protein trajectory.
        """
        if self.set_reference_coordinates(frame_index):
            atom_selection_1 = self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[0]
            )
            atom_selection_2 = self.protein_data.trajectory_data.select_atoms(
                self.atom_selection[1]
            )

            for frame in self.protein_data.frame_indices:
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
    Subclass of ProteinProperty class representing radius of gyration geometric property.
    """

    display_name = "ROG"

    def calculate_property(self):
        """
        Calculate radius of gyration for all trajectory frames for a selection of atoms.
        """
        self.time = []
        for frame in self.protein_data.frame_indices:
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
    Subclass of ProteinProperty representing angles geometric property.

    Parameters
    ----------
    protein_data   : ProteinData
                     An instance of the ProteinData class representing the protein data.
    atom_selection : list
                     List of atom selection with 3 atoms.
    """

    display_name = "A3A"

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
        Calculate distance between 3 given atoms for all trajectory frames.

        Parameters
        ----------
        frame_index : int
                      Reference structure (i.e. frame) from inputed protein trajectory.
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
    Subclass of ProteinProperty class representing dihedral angles geometric property.

    Parameters
    ----------
    protein_data   : ProteinData
                     An instance of the ProteinData class representing the protein data.
    atom_selection : list
                     List of atom selection with 4 atoms.
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
    Subclass of DihedralAngles class representing dihedral angle phi geometric property.
    """

    display_name = "DAPHI4A"

    def calculate_property(self, frame_index=None):
        """
        Calculate dihedral angle phi between 4 given atoms for all trajectory frames.

        Parameters
        ----------
        frame_index : int
                      Reference structure (i.e. frame) from inputed protein trajectory.
        """
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
    Subclass of DihedralAngles class representing dihedral angle psi geometric property.
    """

    display_name = "DAPSI4A"

    def calculate_property(self, frame_index=None):
        """
        Calculate dihedral angle psi between 4 given atoms for all trajectory frames.

        Parameters
        ----------
        frame_index : int
                      Reference structure (i.e. frame) from inputed protein trajectory.
        """
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
