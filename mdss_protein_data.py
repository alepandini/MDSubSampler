import MDAnalysis as mda
import numpy as np


class ProteinData:
    """
    A class used to represent the protein data

    Attributes
    ----------
    trajectory_filename : str
        the path to the trajectory file of the protein
    topology_filename : str
        the path to the topology file of the protein
    config_parameters : str
        configuration parameters of the protein
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
        self.property_dict = {}

    def _read_trajectory(self, trajectory_filename, topology_filename):
        """
        Load trajectory and topology files into Universe to build the object

        Parameters
        ----------------------------
        trajectory_filename: str,
            the path to the trajectory file of the protein.
        topology_filename: str,
            the path to the topology file of the protein

        Returns
        ----------------------------
        Return the number of atoms exist in the object
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
        Read C-Alpha from the first frame of trajectory

        Returns
        ----------------------------
        Return the number of CA atoms from the AtomGroup
        """
        ca_atom_group = self.trajectory_data.select_atoms("name CA")
        return ca_atom_group

    def _frames_of_trajectory(self):
        """
        Method that reads a trajectory and reads the frames that belong to it

        Returns
        ----------------------------
        Return a dictionary that contains frame number, timestep for the frame for
        the whole trajectory
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

        return frames

    def _frame_indices_of_trajectory(self):
        """
        Method that reads a trajectory and reads the frames that belong to it

        Returns
        ----------------------------
        Return a list with the frame indices
        """
        frame_indices = []
        for x in range(len(self.trajectory_data.trajectory)):
            frame_indices.append(x)
        return frame_indices

    def frame_selection_iterator(self, selection_of_frames):
        """
        Method that receives as an input a selection of frames either single frame number
        or sliced selection of frames over a trajectory.

        Returns
        ----------------------------
        A FrameIteratorIndices object with the selected frames. This object has similar
        attributes to a trajectory object and can be used for further analysis.
        """
        trajectory_data = self.trajectory_data.trajectory
        mask = np.array([False for _ in trajectory_data])
        for i in selection_of_frames:
            if isinstance(i, int) or isinstance(i, slice):
                mask[i] = True
            else:
                raise TypeError("Expected int or slice")
        selected_frames = trajectory_data[np.where(mask)[0]]
        return selected_frames

    def frame_selection_indices(self, selection_of_frames):
        """
        Method that receives as an input a selection of frames either single frame number
        or sliced selection of frames over a trajectory.

        Returns
        ----------------------------
        A list with the indices of the selected frames
        """
        trajectory_data = self.trajectory_data.trajectory
        mask = np.array([False for _ in trajectory_data])
        for i in selection_of_frames:
            if isinstance(i, int) or isinstance(i, slice):
                mask[i] = True
            else:
                raise TypeError("Expected int or slice")
        selected_frames = trajectory_data[np.where(mask)[0]]
        indices_of_selected_frames = [ts.frame for ts in selected_frames]
        return indices_of_selected_frames

    def add_property(self, protein_property, property_name):
        self.property_dict[property_name] = protein_property

    def property_data_report(self):
        report_dict = {}
        for k, v in self.property_dict.items():
            report_dict[k] = {
                "min": v.min_value,
                "max": v.max_value,
                "atom_selection": v.atom_selection,
                "name": v.display_name,
            }
        return report_dict
