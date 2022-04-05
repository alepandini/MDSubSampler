import MDAnalysis as mda
import numpy as np

# from pprint import pprint

# import code
# code.interact(local=locals())


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

    def __init__(self, trajectory_filename, topology_filename, config_parameters):

        self.config_par = config_parameters
        self.trajectory_filename = trajectory_filename
        self.topology_filename = topology_filename
        self.trajectory_data = self._read_trajectory(
            self.trajectory_filename, self.topology_filename
        )
        self.n_frames = self.trajectory_data.trajectory.n_frames
        self.ca_atom_group = self._select_CA_atoms()
        self.property_dict = {}
        self.frames = self._frames_of_trajectory(
            self.trajectory_filename, self.topology_filename
        )

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

    def _frames_of_trajectory(self, trajectory, topology):
        """
        Method that reads a trajectory and reads the frames that belong to it

        Returns
        ----------------------------
        Return a dictionary that contains frame number, timestep for the frame for
        the whole trajectory
        """
        trajectory_data = self._read_trajectory(trajectory, topology)
        # frames = {
        #     frame.frame: {"time": frame.time} for frame in trajectory_data.trajectory
        # }
        frames = []
        for x in range(len(trajectory_data.trajectory)):
            _ = trajectory_data.trajectory[x]
            frames.append(
                (
                    x,
                    trajectory_data.trajectory.ts.from_timestep(
                        trajectory_data.trajectory[x]
                    ),
                ),
            )

        return frames

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

    def _select_CA_atoms(self):
        """
        Read C-Alpha from the first frame of trajectory

        Returns
        ----------------------------
        Return the number of CA atoms from the AtomGroup
        """
        ca_atom_group = self.trajectory_data.select_atoms("name CA")
        return ca_atom_group

    def _add_property_dummy(self, protein_property, property_name):
        if property_name in self.property_dict:
            self.property_dict[property_name][sample_label] = protein_property

        else:
            self.property_dict[property_name] = {}
            self.property_dict[property_name][sample_label] = property_property
