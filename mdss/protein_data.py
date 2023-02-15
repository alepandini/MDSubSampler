from datetime import datetime
from re import L
import MDAnalysis as mda
import numpy as np
from mdss.log_setup import log


class ProteinData:
    """
    Represents protein data

    Attributes
    -----------
    trajectory_filename : str
        path to trajectory file
    topology_filename : str
        path to topology file
    config_parameters : str
        protein's configuration parameters
    """

    def __init__(self, trajectory_filename, topology_filename, config_parameters=None):

        self.trajectory_filename = trajectory_filename
        self.topology_filename = topology_filename
        self.trajectory_data = self._read_trajectory(
            self.trajectory_filename, self.topology_filename
        )
        self.ca_atom_group = self._select_CA_atoms()
        self.n_frames = self.trajectory_data.trajectory.n_frames
        self.frames = self._frames_of_trajectory()
        self.frame_indices = self._frame_indices_of_trajectory()
        self.selected_frames = self.frame_selection_iterator([1,5,10,255])
        self.ref_coordinates = self._read_topology(self.topology_filename)
        self.property_dict = {}

    def _read_topology(self, topology_filename):
        top = mda.Universe(topology_filename)
        topology_frame = top.trajectory[0].positions
        return topology_frame

    def _read_trajectory(self, trajectory_filename, topology_filename):
        """
        Load trajectory and topology files into Universe to build the object

        Attributes
        -----------
        trajectory_filename: str,
            path to trajectory file
        topology_filename: str,
            path to topology file

        Returns
        -----------
        Number of atoms exist in object
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
        Reads C-Alpha from the first frame of trajectory

        Returns
        -----------
        Number of CA atoms from AtomGroup
        """
        ca_atom_group = self.trajectory_data.select_atoms("name CA")
        return ca_atom_group

    def _frames_of_trajectory(self):
        """
        Creates a dictionary with frame number and timestep for full trajectory
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
        Creates list with frames' indices for full trajectory
        """
        frame_indices = []
        for x in range(len(self.trajectory_data.trajectory)):
            frame_indices.append(x)
        return frame_indices

    def frame_selection_iterator(self, selection_of_frames):
        """
        Attributes
        -----------
        selection_of_frames: int,
            single frame or slice of frames from trajectory

        Returns
        -----------
        FrameIteratorIndices object with the selected frames: This object has similar
        attributes to a trajectory object and can be used for further analysis.
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
        Attributes
        -----------
        selection_of_frames: int,
            single frame or slice of frames from trajectory

        Returns
        -----------
        List with indices of selected frames
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
        log.info(
            "{:15s} Indices of selected frames: {}".format(
                "OUTPUT", indices_of_selected_frames
            )
        )
        return indices_of_selected_frames

    def write_xtc_file(self):
        protein = self.trajectory_data.select_atoms("protein")
        selection = self.selected_frames
        with mda.Writer("test_protein.xtc", protein.n_atoms) as W:
	        for ts in selection:
		        W.write(protein)

    # def cast_output_ML():
    #     # selection_of_frames = 
    #     # selected_frames_trajectory = frame_selection_iterator(self, selection_of_frames)
    #     pass

    def add_property(self, protein_property, property_name):
        """
        Gets key from property dictionary
        """
        timestamp = str(datetime.now().timestamp())
        key = "{}_{}".format(property_name, timestamp)
        self.property_dict[key] = protein_property
        return key

    def property_data_report(self):
        """
        Creates a report with key information and statistics for property
        """
        report_dict = {}
        for k, v in self.property_dict.items():
            report_dict[k] = {
                "min": v.min_value,
                "max": v.max_value,
                "atom_selection": v.atom_selection,
                "name": v.display_name,
            }
        return report_dict
