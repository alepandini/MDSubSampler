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
from sklearn.model_selection import train_test_split
from datetime import datetime
from copy import deepcopy
from re import L
import MDAnalysis as mda
import numpy as np
from mdss.log_setup import log
import numpy as np
import json


class ProteinData:
    """
    Class representing protein data.

    Attributes
    -----------
    trajectory_filename : str
        Path to trajectory file.
    topology_filename : str
        Path to topology file.
    config_parameters : str, optional
        Protein's configuration parameters. Default is None.
    """

    def __init__(
        self,
        trajectory_filename,
        topology_filename,
        config_parameters=None,
    ):
        """
        Initialise ProteinData object.

        Parameters
        -----------
        trajectory_filename : str
            Path to trajectory file.
        topology_filename : str
            Path to topology file.
        config_parameters : str, optional
            Protein's configuration parameters. Default is None.
        """

        self.trajectory_filename = trajectory_filename
        self.topology_filename = topology_filename
        self.trajectory_data = self._read_trajectory(
            self.trajectory_filename, self.topology_filename
        )
        self.topology_data = self._read_topology(self.topology_filename)
        self.ca_atom_group = self._select_CA_atoms()
        self.n_frames = self.trajectory_data.trajectory.n_frames
        self.frames = self._frames_of_trajectory()
        self.frame_indices = self._frame_indices_of_trajectory()
        self.ref_coordinates = self.topology_data.trajectory[0].positions
        self.property_dict = {}

    def _read_topology(self, topology_filename):
        """
        Load topology file.

        Parameters
        ----------
        topology_filename : str
            Path to topology file.

        Returns
        -------
        mda.Universe
            An instance of the MDAnalysis Universe representing the loaded topology data.
        """
        top_data = mda.Universe(topology_filename)
        return top_data

    def _read_trajectory(self, trajectory_filename, topology_filename):
        """
        Load trajectory and topology files into Universe to build the object.

        Parameters
        -----------
        trajectory_filename : str
            Path to trajectory file.
        topology_filename : str
            Path to topology file.

        Returns
        -----------
        mda.Universe
            An instance of the MDAnalysis Universe representing the loaded trajectory.
        """
        trajectory_data = mda.Universe(
            topology_filename,
            trajectory_filename,
            permissive=False,
            topology_format="GRO",
        )
        return trajectory_data

    def _select_CA_atoms(self):
        """
        Select C-alpha atoms from the first frame of the trajectory.

        Returns
        -------
        MDAnalysis.core.groups.AtomGroup
            An AtomGroup containing the C-alpha atoms from the first frame of the trajectory.
        """
        ca_atom_group = self.trajectory_data.select_atoms("name CA")
        return ca_atom_group

    def _frames_of_trajectory(self):
        """
        Generate a dictionary with frame numbers and timesteps for a protein trajectory.

        Returns
        -------
        list of tuples
            A list of tuples containing the frame number (index) and corresponding timestep for each frame.
        """
        frames = []
        for x in range(len(self.trajectory_data.trajectory)):
            _ = self.trajectory_data.trajectory[x]
            frames.append(
                (
                    x,
                    self.trajectory_data.trajectory.ts.from_timestep(
                        self.trajectory_data.trajectory[x]
                    ),
                ),
            )

        log.info("{:15s} Number of frames: {}".format("INPUT", len(frames)))
        return frames

    def _frame_indices_of_trajectory(self):
        """
        Generate a list of frame indices for a protein trajectory.

        Returns
        -------
        list
            A list of integers representing the frame indices from the protein trajectory.

        """
        frame_indices = []
        for x in range(len(self.trajectory_data.trajectory)):
            frame_indices.append(x)
        return frame_indices

    def frame_selection_iterator(self, selection_of_frames):
        """
        Create a new object with similar attributes to a trajectory object from a specific selection of frames.

        Parameters
        ----------
        selection_of_frames : int or slice
            Single frame or slice of frames from the trajectory to select.

        Returns
        -------
        FrameIteratorIndices
            An instance of the MDAnalysis.coordinates.base.FrameIteratorIndices.
            It is iterable over the frames of a trajectory.

        Raises
        ------
        TypeError
            If the `selection_of_frames` parameter is neither an integer nor a slice.

        Notes
        -----
        The method creates a boolean mask array to indicate the selected frames.
        If an integer or slice is provided, the corresponding indices in the mask are set to True.
        The selected frames are extracted from the trajectory data using the mask.
        """
        trajectory_data = self.trajectory_data.trajectory
        mask = np.array([False for _ in trajectory_data])
        for i in selection_of_frames:
            if isinstance(i, int) or isinstance(i, slice):
                mask[i] = True
            else:
                log.error(
                    "{:15s} Expected int or slice in frame_selection_iterator()".format(
                        "STEPS"
                    )
                )
                raise TypeError("Expected int or slice")
        selected_frames = trajectory_data[np.where(mask)[0]]
        log.info(
            "{:15s} Number of selected frames: {}".format(
                "OUTPUT", len(selected_frames)
            )
        )
        return selected_frames

    def frame_selection_indices(self, selection_of_frames):
        """
        Generate a list with only selected frames from a protein trajectory

        Parameters
        -----------
        selection_of_frames : int or slice
            Single frame or slice of frames from the trajectory to select.

        Returns
        -------
        List
            Contains indices of selected frames.

        Raises
        ------
        TypeError
            If the `selection_of_frames` parameter is neither an integer nor a slice.
        """
        trajectory_data = self.trajectory_data.trajectory
        mask = np.array([False for _ in trajectory_data])
        for i in selection_of_frames:
            if isinstance(i, int) or isinstance(i, slice):
                mask[i] = True
            else:
                log.error(
                    "{:15s} Expected int or slice in frame_selection_iterator()".format(
                        "STEPS"
                    )
                )
                raise TypeError("Expected int or slice")
        selected_frames = trajectory_data[np.where(mask)[0]]
        indices_of_selected_frames = [ts.frame for ts in selected_frames]
        return indices_of_selected_frames

    def write_xtc_file(self, outfilepath, selected_frames):
        """
        Generate an xtc file containing only selected frames from a protein trajectory.

        Parameters
        -----------
        outfilepath : str
            Path where output file is saved.
        selected_frames : int ot list,
            Single frame or list of frames from trajectory.
        """
        protein = self.trajectory_data.select_atoms("protein")
        with mda.Writer(outfilepath, protein.n_atoms) as W:
            for t_idx in selected_frames:
                self.trajectory_data.trajectory[t_idx]
                W.write(protein)

    def cast_output_traj_to_numpy(self, outfilepath, subsampled_traj, unit="nanometer"):
        """
        Casts an XTC file into a NumPy array for user readability.

        Parameters
        -----------
        outfilepath : str
            Path where output file is saved.
        subsampled_traj : MDAnalysis.coordinates.XTC.XTCReader
            XTC trajectory file.
        unit : str, optional
            Unit for coordinates valuess.

        Returns
        -------
        numpy.ndarray
            NumPy array containing the coordinates of the subsampled trajectory.

        """
        coordinates_numpy = []
        for ts in subsampled_traj:
            coordinates_numpy.append(deepcopy(ts.positions))
        coordinates_numpy = np.array(coordinates_numpy)
        if unit == "nanometer":
            coordinates_numpy = coordinates_numpy / 10
        np.save(outfilepath, coordinates_numpy)
        return coordinates_numpy

    def convert_numpy_to_2D(self, infilepath, outfilepath):
        """
        Convert a 3D numpy array to a 2D numpy array and save it to a file.

        Parameters
        ----------
        infilepath : numpy.ndarray
            The input 3D numpy array to be converted.
        outfilepath : str
            The path where the output file will be saved.

        Returns
        -------
        numpy.ndarray
            The converted 2D numpy array.
        """
        (x, y, z) = infilepath.shape
        outfile = np.reshape(infilepath, (x, y * z))
        np.save(outfilepath, outfile)
        return outfile

    def ML_input_prep(self, infilepath, outfilepath_training, outfilepath_testing):
        """
        Prepares input data for machine learning by splitting the input file into training and testing data.

        Parameters
        ----------
        infilepath : str
            Path to the input file containing the data to be split.
        outfilepath_training : str
            Path where the training data file will be saved.
        outfilepath_testing : str
            Path where the testing data file will be saved.
        """
        training_data, testing_data = train_test_split(
            infilepath, test_size=0.3, random_state=25
        )
        np.save(outfilepath_training, training_data)
        np.save(outfilepath_testing, testing_data)

    def add_property(self, protein_property, property_name):
        """
        Add a protein property to the dictionary.

        Parameters
        ----------
        protein_property : ProteinProperty
            An object of the ProteinProperty class that represents the protein property.
        property_name : str
            The name of the property to be added.

        Returns
        -------
        str
            A string containing the property name and timestamp.
        """
        timestamp = str(datetime.now().timestamp())
        key = "{}_{}".format(property_name, timestamp)
        self.property_dict[key] = protein_property
        return key

    def property_data_report(self, outfilepath):
        """
        Create a JSON report with key information and statistics for the property.

        Parameters
        ----------
        outfilepath : str
            Path to the output file where the JSON report will be saved.

        Returns
        -------
        dict
            A dictionary containing the report information.
        """
        report_dict = {}
        for k, v in self.property_dict.items():
            report_dict[k] = {
                "min": round(v.min_value, 3),
                "max": round(v.max_value, 3),
                "atom_selection": v.atom_selection,
                "property_name": v.display_name,
                "dissimilarity": round(v.ref_dissimilarity, 5),
                "traj_size": len(v.frame_indices),
            }
            if hasattr(v, "ref_property"):
                report_dict[k]["sample_percent"] = (
                    100 * len(v.frame_indices) / len(v.ref_property.frame_indices)
                )
        with open(outfilepath, "w") as f:
            json.dump(report_dict, f, indent=2)

        return report_dict
