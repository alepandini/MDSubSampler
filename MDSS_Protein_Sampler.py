
import MDAnalysis as mda
import numpy as np
import pandas as pd
import random
import dictances
import pingouin as pg

from MDAnalysis.analysis import rms
from Protein_Sampler import discretize_to_dict, replace_zero


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

    def output_trj_summary(self):
        print("----------------------")
        print("TRAJECTORY INFORMATION")
        print("----------------------")
        # Printing total no of frames , atoms and C-Alpha in trajectroy
        print(
            "n frames = {0}\nn atoms = {1}\nn CA atoms = {2}".format(
                self.n_frames,
                self.trajectory_data.trajectory.n_atoms,
                self.ca_atom_group.n_atoms,
            )
        )

    def statistical_analysis(self, stat_vector, stat_sampled_vector, sample_label):

        print("----------------------")
        print("STATISTICAL ANALYSIS")
        print("----------------------")
        vector_mean = np.mean(np.array([stat_vector]), axis=0)
        vector_variance = np.var(stat_vector)
        print("mean for {0} vector::{1}".format(sample_label, vector_mean))
        print("variance for {0} vector::{1}".format(sample_label, vector_variance))
        print("----------------------")
        print("Kolmogorov-Smirnov test")
        print("----------------------")
        # perform Kolmogorov-Smirnov test on two vectors
        KSstat, pvalue = ks_2samp(stat_vector, stat_sampled_vector)
        print("stat=%.3f, p=%.3f" % (KSstat, pvalue))
        # interpret the Kolmogorov-Smirnov test results
        if pvalue > 0.05:
            print("Probably the same distribution")
        else:
            print("Probably different distributions")
        print("----------------------")
        print("Student T-test")
        print("----------------------")
        # perform t-test on two vectors
        print(pg.ttest(x=stat_vector, y=stat_sampled_vector, correction=True).round(2))
        stat, p1 = ttest_ind(stat_vector, stat_sampled_vector)
        print("stat1=%.3f, p1=%.3f" % (stat, p1))
        # interpret the t-test results
        if p1 > 0.05:
            print("Probably the same distribution")
        else:
            print("Probably different distributions")

    def add_property(self, property_vector, property_name, sample_label):
        if property_name in self.property_dict:
            self.property_dict[property_name][sample_label] = property_vector

        else:
            self.property_dict[property_name] = {}
            self.property_dict[property_name][sample_label] = property_vector


# Added the vector argument so I can calculate min, max, avg and do discretize_to_dict
class ProteinProperty:
    def __init__(self, protein_data, vector, frame_list, atom_selection="name CA"):
        self.protein_data = protein_data
        self.atom_selection = atom_selection
        self.property_vector = []

    def discretize_vector(self):
        self.property_vector_discretized = discretize_to_dict(self.property_vector, self.min_value, self.max_value)

    def _property_statistics(self):
        self.min_value = np.min(self.property_vector)
        self.max_value = np.max(self.property_vector)
        self.avg_value = np.average(self.property_vector) 

    def set_reference_coordinates(self):
        self.protein_data.trajectory_data.trajectory[0]  # setting us on the first frame
        self.ref_coordinates = self.protein_data.trajectory_data.select_atoms(
            self.atom_selection
        ).positions.copy()  # extracting a copy of the coordinates of the first frame only for a selection of atoms

class RMSDProperty(ProteinProperty):
    def __init__(self, protein_data, frame_list, atom_selection="name CA"):

        super().__init__(protein_data, frame_list, atom_selection)

        self.set_reference_coordinates()
        
        for frame in frame_list:
            # go through the trajectory and for each frame I compare with my reference frame
            self.protein_data.trajectory_data.trajectory[frame]
            self.property_vector.append(                                                              
                rms.rmsd(
                    self.protein_data.trajectory_data.select_atoms(atom_selection).positions,
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
            self.property.append(
                protein_data.trajectory_data.select_atoms(
                    atom_selection
                ).radius_of_gyration()
            )


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
        return dictances.bhattacharyya(self.property_1.property_vector_discretized, self.property_2.property_vector_discretized)


class KLDiverDistance(Distance):
    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        return dictances.kullback_leibler(self.property_1.property_vector_discretized, self.property_2.property_vector_discretized)


class PearsonDictDistance(Distance):
    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_distance(self):
        return dictances.pearson(self.property_1.property_vector_discretized, self.property_2.property_vector_discretized)


### # replace 0 with small number and then rescale values to make the sum equal to 1
### if clean:
###     self.prop_vector1 = replace_zero(self.freq_prop_vector1)
###     self.prop_vector1 = replace_zero(self.freq_prop_vector1)

### print("size: {0:6d} distance: {1:.3f}".format(sample_size, distance))
### 
### print("-------------------------")
### 
### for sample_size in [10, 50, 100, 200, 400, 500]:
###     # for sample_size in [800,4000,8000,16000,32000,40000]:
### 
###     sub_prop_vector = random.sample(Prop_vector, sample_size)
### 
###     # discretisation of the subsampled vector
###     freq_sub_prop_vector = discretize_to_dict(
### 	sub_prop_vector, min_value, max_value
###     )
### 
###     if clean:
### 	freq_sub_prop_vector = replace_zero(freq_sub_prop_vector)
###     # calculate the kl divergence
###     pq_distance = self.calculate_distance(
### 	freq_sub_prop_vector, freq_sub_prop_vector
###     )
### 
###     print("size: {0:6d} distance: {1:.3f}".format(sample_size, pq_distance))
