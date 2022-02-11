from dataclasses import replace
import MDAnalysis as mda
import numpy as np
import pandas as pd
import random
import dictances

from pprint import pprint

# import code
# code.interact(local=locals())

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.dihedrals import Ramachandran


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


class ProteinProperty:
    """
    A class used to calculate a property for the protein based on the frame selection
    from the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    vector:
        Simple vector used to calculate statistics
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = None

    def __init__(self, protein_data, vector, frame_list, atom_selection="name CA"):
        self.protein_data = protein_data
        self.atom_selection = atom_selection
        self.property_vector = []

    def _add_reference_to_protein_data(self):
        """
        Method that links the ProteinProperty and ProteinData classes
        """
        self.protein_data._add_property_dummy(self, self.property_name)

    def discretize_vector(self, min_value=None, max_value=None):
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
        bin_size = (max_value - min_value) / 100.0
        bin_vector = np.arange(min_value, max_value, bin_size)
        counts, bins = np.histogram(self.property_vector, bins=bin_vector)
        self.property_vector_discretized = dict(
            zip(bins, (counts / len(self.property_vector)))
        )
        return self.property_vector_discretized

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
        self.protein_data.trajectory_data.trajectory[0]
        self.ref_coordinates = self.protein_data.trajectory_data.select_atoms(
            self.atom_selection
        ).positions.copy()

    def write_property_vector(self, outfilename):
        """
        Method that saves the vector with the calculations of a specific property for
        a protein in a file
        """
        with open(outfilename, "w") as f:
            for i, value in enumerate(self.property_vector):
                f.write("{} {}\n".format(i, value))

    def write_property_discretised_vector(self, outfilename):
        """
        Method that saves the discretised vector with the calculations of a specific
        property for a protein in a file
        """
        discr_vector = self.discretize_vector()
        with open(outfilename, "w") as f:
            for (key, value) in discr_vector.items():
                f.write("{} {}\n".format(key, value))


class RMSDProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the RMSD value for each frame in the
    protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = "RMSD"

    def __init__(self, protein_data, frame_list, atom_selection="name CA"):
        super().__init__(protein_data, frame_list, atom_selection)
        self.set_reference_coordinates()

        for frame in frame_list:
            """
            Go through the trajectory and for each frame I compare with my reference frame
            """
            self.protein_data.trajectory_data.trajectory[frame]
            self.property_vector.append(
                rms.rmsd(
                    self.protein_data.trajectory_data.select_atoms(
                        atom_selection
                    ).positions,
                    self.ref_coordinates,
                )
            )

        self._property_statistics()
        self.discretize_vector()


class DistanceProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the distance value between two
    atoms or two group of atoms for each frame in the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    # def __init__(self, protein_data, frame_list, atom_selection="name CA"):

    # super().__init__(protein_data, frame_list, atom_selection)
    pass


class RMSFProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the RMSF value for the particles
    in the atomgroup in each frame in the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = "RMSF"

    def __init__(self, protein_data, frame_list, atom_selection="name CA"):
        super().__init__(protein_data, frame_list, atom_selection)
        self.set_reference_coordinates()

        u = self.protein_data.trajectory_data
        protein = u.select_atoms("protein")
        # probably changes the u so maybe I need to make a copy of it
        prealigner = align.AlignTraj(u, u, select=atom_selection, in_memory=True).run()
        reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
        reference = mda.Merge(protein).load_new(
            reference_coordinates[:, None, :], order="afc"
        )
        aligner = align.AlignTraj(
            u, reference, select="protein and name CA", in_memory=True
        ).run()
        ca = protein.select_atoms(atom_selection)
        print(ca)
        rmsfer = rms.RMSF(ca, verbose=True).run()
        self.property_vector.append(rmsfer.results.rmsf)

        self._property_statistics()
        self.discretize_vector()


class RadiusOfGyrationProperty(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the Radius of Gyration value for each frame
    in the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = "rog"

    def __init__(self, protein_data, frame_list, atom_selection="name CA"):

        super().__init__(protein_data, frame_list, atom_selection)

        self.time = []
        for frame in frame_list:
            protein_data.trajectory_data.trajectory[frame]

            self.time.append(protein_data.trajectory_data.trajectory.time)
            self.property_vector.append(
                protein_data.trajectory_data.select_atoms(
                    atom_selection
                ).radius_of_gyration()
            )

        self._property_statistics()
        self.discretize_vector()


class DihedralAngles(ProteinProperty):
    """
    A Subclass of ProteinProperty class used to calculate the angles for selected set of
    atoms through the protein trajectory

    Attributes
    ----------
    protein_data : ProteinData object
        Contains the trajectory and topology data for the protein
    frame_list: list
        List that contains all the frames from a given protein trajectory
    atom_selection: str
        Choice of atoms for calculation of a property on this selection of atoms
    """

    display_name = "dihedral"

    def __init__(self, protein_data, frame_list, atom_selection="name CA"):
        super().__init__(protein_data, frame_list, atom_selection)
        self.set_reference_coordinates()

        for frame in frame_list:
            self.protein_data.trajectory_data.trajectory[frame]
            R = Ramachandran(atom_selection).run()
            self.property_vector.append(R.Ramachandran.angles)

        self._property_statistics()
        self.discretize_vector()


class ProteinSampler:
    """
    A class used to create a sample of a protein trajectory

    Attributes
    ----------
    frame_list: list
        List that contains all the frames from a given protein trajectory
    """

    def __init__(self, frame_list):
        self.frame_list = frame_list
        self.sampled_frame_list = None

    def sample(self, size):
        self.sampled_frame_list.sort()


class RandomSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Random Sampling

    Attributes
    ----------
    frame_list: list
        List that contains all the frames from a given protein trajectory
    seed: int
        Number that initialise a random-number generator
    """

    def __init__(self, frame_list, seed_number=1999):
        random.seed(seed_number)
        super().__init__(frame_list)

    def sample(self, size):
        """
        Method that generates a random sample of a list

        Attributes
        ----------
        size: int
            The sample size

        Returns
        ----------
        return a sample of the frame list with the desired size
        """
        self.sampled_frame_list = random.sample(self.frame_list, size)
        super().sample(size)
        return self.sampled_frame_list


class StratifiedSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Stratified Sampling
    """

    def __init__(self, frame_list, layers):
        self.layers = layers
        super().__init__(frame_list)

    def strata_sample_size(self, size, population, layer_size):
        """
        Method that calculates the sample size of the strata (ie homogeneous groups)

        Attributes
        ----------
        size: int
            Whole sample size
        population: int
            Whole population size (sum of each layer - ie strata - size)
        layer_size: int
            Size for current layer

        Returns
        ----------
        return the rounded sample size of the strata

        """
        cur_size = (size / population) * layer_size
        return round(cur_size)

    def sample(self, size):
        """
        Method that does the stratified sampling

        Attributes
        ----------
        layers: vector
            This is a 2D vector that consists of multiple layers. Each layer is a set of labels
            for the frames according to the strata
        size: int
            Whole sample size

        Returns
        ----------
        return a list of samples

        """
        population = sum(len(layer) for layer in self.layers)
        samples = []

        for layer in self.layers:
            layer_size = len(layer)
            current_layer_sample_size = self.strata_sample_size(
                size, population, layer_size
            )
            current_sample = random.sample(layer, current_layer_sample_size)
            samples.extend(current_sample)

        return samples


class BootstrappingSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Bootstrapping Sampling

    Attributes
    ----------
    frame_list: list
        List that contains all the frames from a given protein trajectory
    """

    def sample(self, size, number_of_iterations):
        """
        Method that does the bootstrapping sampling

        Attributes
        ----------
        size: int
            This is the desired size of the sample each time we iterate
        number_of_iterations: int
            This is the number of times the random sampling method is performed

        Returns
        ----------
        return a list of samples

        """
        samples = []
        for i in range(number_of_iterations):
            current_sample = random.sample(self.frame_list, size)
            current_sample_mean = self.find_nearest(
                current_sample, np.mean(current_sample)
            )
            samples.append(current_sample_mean)

        return samples

    def find_nearest(self, array, value):
        """
        Method that finds the closest number to the mean value

        Attributes
        ----------
        array: list
            This is a list with the data
        value: float
            This is the number of times the random sampling method is performed

        Returns
        ----------
        return an integer with the frame number that is closest to the mean number

        """
        idx = (np.abs(array - value)).argmin()
        return array[idx]


# vector_1 for full protein and vector_2 for sample
class Distance:
    """
    A class used to calculate the distance in terms of property between
    a full protein trajectory and a sample of it to identify the difference

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = None

    def __init__(self, property_1, property_2, clean=False):
        self.property_1 = property_1
        self.property_2 = property_2
        self.distance = self.calculate_distance()

    def calculate_distance(self):
        """
        Method that calculates the difference between the average values of the
        two calculated property vectors.
        """
        return self.property_1.avg_value - self.property_2.avg_value


class BhattaDistance(Distance):
    """
    A Subclass of the Distance class that calculates the Bhattacharya distance between
    two property vectors

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = "bhatta"

    def __init__(self, property_1, property_2):
        self.min_value = min(
            min(property_1.property_vector), min(property_2.property_vector)
        )
        self.max_value = max(
            max(property_1.property_vector), max(property_2.property_vector)
        )
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        """
        Method that returns the Bhatta distance between two vectors
        """
        property_1_discretized = self.property_1.discretize_vector(
            min_value=self.min_value, max_value=self.max_value
        )
        property_2_discretized = self.property_2.discretize_vector(
            min_value=self.min_value, max_value=self.max_value
        )

        return dictances.bhattacharyya(property_1_discretized, property_2_discretized)


class KLDiverDistance(Distance):
    """
    A Subclass of the Distance class that calculates the Kullback-Leibler divergence between
    two property vectors

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = "kl"

    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        """
        Method that returns the KL distance between two vectors
        """
        return dictances.kullback_leibler(
            self.property_1.property_vector_discretized,
            self.property_2.property_vector_discretized,
        )


class PearsonDictDistance(Distance):
    """
    A Subclass of the Distance class that calculates the Pearson distance between
    two property vectors

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = "pearson"

    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        """
        Method that returns the pearson distance between two vectors
        """
        return dictances.pearson(
            self.property_1.property_vector_discretized,
            self.property_2.property_vector_discretized,
        )
