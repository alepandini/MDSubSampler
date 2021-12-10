import MDAnalysis as mda
import numpy as np
import pandas as pd
import random
import dictances

from MDAnalysis.analysis import rms


class ProteinData:
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
        # Load trajectory file into Universe
        trajectory_data = mda.Universe(
            topology_filename,
            trajectory_filename,
            permissive=False,
            topology_format="GRO",
        )
        return trajectory_data

    def _select_CA_atoms(self):
        # Read C-Alpha from the first frame of trajectory
        ca_atom_group = self.trajectory_data.select_atoms("name CA")
        return ca_atom_group

    def _add_property_dummy(self, protein_property, property_name):
        if property_name in self.property_dict:
            self.property_dict[property_name][sample_label] = protein_property

        else:
            self.property_dict[property_name] = {}
            self.property_dict[property_name][sample_label] = property_property


# Added the vector argument so I can calculate min, max, avg and do discretize_to_dict
class ProteinProperty:
    def __init__(self, protein_data, vector, frame_list, atom_selection="name CA"):
        self.protein_data = protein_data
        self.atom_selection = atom_selection
        self.property_vector = []

    def _add_reference_to_protein_data(self):
        self.protein_data._add_property_dummy(self, self.property_name)

    def discretize_vector(self):
        bin_size = (self.max_value - self.min_value) / 100.0
        bin_vector = np.arange(self.min_value, self.max_value, bin_size)
        counts, bins = np.histogram(self.property_vector, bins=bin_vector)
        self.property_vector_discretized = dict(
            zip(bins, counts / len(self.property_vector))
        )

    def _property_statistics(self):
        self.min_value = np.min(self.property_vector)
        self.max_value = np.max(self.property_vector)
        self.avg_value = np.average(self.property_vector)

    def set_reference_coordinates(self):
        self.protein_data.trajectory_data.trajectory[0]  # setting us on the first frame
        self.ref_coordinates = self.protein_data.trajectory_data.select_atoms(
            self.atom_selection
        ).positions.copy()  # extracting a copy of the coordinates of the first frame only for a selection of atoms

    def write_property_vector(self, outfilename):
        fout = open(outfilename, "w")
        for value in self.property_vector:
            fout.write("{0}\n".format(value))
        fout.close()


class RMSDProperty(ProteinProperty):
    def __init__(self, protein_data, frame_list, atom_selection="name CA"):

        super().__init__(protein_data, frame_list, atom_selection)

        self.set_reference_coordinates()

        for frame in frame_list:
            # go through the trajectory and for each frame I compare with my reference frame
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


class RadiusOfGyrationProperty(ProteinProperty):
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


class ProteinSampler:
    def __init__(self, frame_list):
        self.frame_list = frame_list
        self.sampled_frame_list = None

    def sample(self, size):
        self.sampled_frame_list.sort()


class RandomSampler(ProteinSampler):
    def __init__(self, frame_list, seed_number=1999):
        random.seed(seed_number)
        super().__init__(frame_list)

    def sample(self, size):
        self.sampled_frame_list = random.sample(self.frame_list, size)
        super().sample(size)
        return self.sampled_frame_list


class StratifiedSampler(ProteinSampler):

    # size = whole sample size
    # population = whole population size (sum of each layer size)
    # layer_size = size for current layer
    def strata_sample_size(size, population, layer_size):
        cur_size = (size / population) * layer_size

        return round(cur_size)

    # strata_vector = 2D vector
    # the strata vector consists of multiple layers
    # each layer is a set of labels for the frames according to the strata
    # size = whole sample size
    def stratified_sampling(layers, size):
        population = sum(len(layer) for layer in layers)
        samples = []

        for layer in layers:
            layer_size = len(layer)
            print(f"layer size: {layer_size}")
            current_layer_sample_size = self.strata_sample_size(
                size, population, layer_size
            )
            print(f"layer sample size: {current_layer_sample_size }")
            print(f"proportion: {current_layer_sample_size / layer_size}")
            current_sample = random.sample(layer, current_layer_sample_size)
            print(f"sample: {current_sample}")
            samples.extend(current_sample)

        return samples


# vector_1 for full protein and vector_2 for sample
class Distance:
    def __init__(self, property_1, property_2, clean=False):
        self.property_1 = property_1
        self.property_2 = property_2
        self.distance = self.calculate_distance()

    def calculate_distance(self):
        return self.property_1.avg_value - self.property_2.avg_value


class BhattaDistance(Distance):
    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        return dictances.bhattacharyya(
            self.property_1.property_vector_discretized,
            self.property_2.property_vector_discretized,
        )


class KLDiverDistance(Distance):
    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        return dictances.kullback_leibler(
            self.property_1.property_vector_discretized,
            self.property_2.property_vector_discretized,
        )


class PearsonDictDistance(Distance):
    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        return dictances.pearson(
            self.property_1.property_vector_discretized,
            self.property_2.property_vector_discretized,
        )
