import numpy as np
import MDAnalysis as mda

from mdss_protein_data import ProteinData


class ProteinProperty:
    """
    A class used to calculate a property for the protein based on the frame selection
    from the protein trajectory

    Attributes
    ----------
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class

    vector:
        Simple vector used to calculate statistics

    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = None

    def __init__(self, protein_data, atom_selection="name CA"):
        if not isinstance(protein_data, ProteinData):
            print("Warning: A instance of ProteinData is required. protein_data attribute set to None")
            protein_data = None
        self.protein_data = protein_data
        self.atom_selection = atom_selection
        self.ref_coordinates = []
        self.property_vector = []
        self.property_vector_discretized = {}
        self.frame_indices = []
        self._add_reference_to_protein_data()

    @classmethod
    def from_xvg(cls, xvg_filepath):
        instance = cls(protein_data=None, atom_selection=[None, None])
        _frames, distance_values = np.loadtxt(xvg_filepath, unpack=True)
        instance.property_vector = distance_values
        instance._property_statistics()
        instance.discretize_vector()
        return instance

    def _add_reference_to_protein_data(self):
        """
        Method that links the ProteinProperty and ProteinData classes
        """
        if self.protein_data is not None:
            self.protein_data.add_property(self, self.display_name)

    def discretize_vector(self, min_value=None, max_value=None, n_bins=100):
        """
        Method that discretises a vector used for Bhatta and KL distances

        Returns
        ----------------------------
        return the discretised vector
        """
        if min_value is None:
            min_value = self.min_value
        if max_value is None:
            max_value = self.max_value
        bin_size = (max_value - min_value) / n_bins
        bin_vector = np.arange(min_value, max_value, bin_size)
        counts, bins = np.histogram(self.property_vector, bins=bin_vector)
        self.property_vector_discretized = dict(
            zip(bins, (counts / len(self.property_vector)))
        )

    def _property_statistics(self):
        """
        Method that calculates the minimum, maximum and average values of a vector
        """
        self.min_value = np.min(self.property_vector)
        self.max_value = np.max(self.property_vector)
        self.avg_value = np.average(self.property_vector)

    def set_reference_coordinates(self):
        """
        Method that sets us on the first frame and extracts a copy of the coordinates
        of the first frame only for a selection of atoms
        """
        if self.protein_data is not None:
            self.protein_data.trajectory_data.trajectory[0]
            if isinstance(self.atom_selection, list):
                self.ref_coordinates = [
                    self.protein_data.trajectory_data.select_atoms(
                        selection
                    ).positions.copy()
                    for selection in self.atom_selection
                ]
            else:
                self.ref_coordinates = self.protein_data.trajectory_data.select_atoms(
                    self.atom_selection
                ).positions.copy()
            return True
        else:
            print('Warning: property is not associated to a protein data object.')
            return False

    def write_property_vector(self, outfilepath):
        """
        Method that saves the vector with the calculations of a specific property for
        a protein in a file
        """
        with open(outfilepath, "w") as f:
            for i, value in enumerate(self.property_vector):
                f.write("{} {}\n".format(i, value))

    def write_property_discretised_vector(self, outfilepath):
        """
        Method that saves the discretised vector with the calculations of a specific
        property for a protein in a file
        """
        discr_vector = self.discretize_vector()
        with open(outfilepath, "w") as f:
            for (key, value) in discr_vector.items():
                f.write("{} {}\n".format(key, value))


class SampledProperty(ProteinProperty):

    display_name = "Sampled_Property"

    def __init__(self, property_vector, indices):
        self.property_vector = property_vector
        self.indices = indices
        self._property_statistics()

    def calculate_property(self):
        pass

    def write_property_vector_sample(self, outfilepath):
        """
        Method that saves the vector with the calculations of a specific property for
        a protein in a file
        """
        with open(outfilepath, "w") as f:
            for i, value in zip(self.indices, self.property_vector):
                f.write("{} {}\n".format(i, value))
