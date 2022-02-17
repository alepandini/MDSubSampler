import MDAnalysis as mda
from pprint import pprint

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

    def _add_property_dummy(self, protein_property, property_name):
        if property_name in self.property_dict:
            self.property_dict[property_name][sample_label] = protein_property

        else:
            self.property_dict[property_name] = {}
            self.property_dict[property_name][sample_label] = property_property

    def selection_of_frames(self, trajectory, topology):
        """
        Frame selection from the trajectory

        Returns
        ----------------------------
        Return a dictionary that contains frame number, timestep for the frame
        """
        trajectory_data = self._read_trajectory(trajectory, topology)
        frames = {
            frame.frame: {"time": frame.time} for frame in trajectory_data.trajectory
        }
        return frames
