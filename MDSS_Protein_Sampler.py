import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca
import numpy.linalg
import numpy as np
import pandas as pd
import os
import configparser
import random
import dictances
import itertools
import seaborn as sns
import scipy
import pingouin as pg
import csv
import matplotlib.pyplot as plt


from sklearn import preprocessing
from scipy import stats
from pandas.plotting import scatter_matrix
from matplotlib import pyplot
from MDAnalysis.analysis import rms
from MDAnalysis.coordinates.memory import MemoryReader
from dictances import bhattacharyya, bhattacharyya_coefficient, kullback_leibler
from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction, analysis_class
from scipy.stats import mannwhitneyu, shapiro, ks_2samp, ttest_ind
from scipy.signal import savgol_filter


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
        # self.masses = self.ca_atom_group.masses
        # self.total_mass = np.sum(self.masses)
        self.property_dict = {}

    def _read_trajectory(self, trajectory_filename, topology_filename):
        # Load trajectory file into Universe
        trajectory_data = mda.Universe(
            topology_filename,
            trajectory_filename,
            permissive=False,
            topology_format="PDB",
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


class ProteinProperty:
    """TODO.

    Attributes:
        property TODO
    """

    pass


class RMSDProperty(ProteinProperty):
    def __init__(self, protein_data, frame_list, atom_selection="name CA"):
        protein_data.trajectory_data.trajectory[0]
        ref_coordinates = protein_data.trajectory_data.select_atoms(
            atom_selection
        ).positions.copy()

        self.property = []
        for frame in frame_list:
            protein_data.trajectory_data.trajectory[frame]
            self.property.append(
                rms.rmsd(
                    protein_data.trajectory_data.select_atoms(atom_selection).positions,
                    ref_coordinates,
                )
            )


class RadiusOfGyrationProperty(ProteinProperty):
    def __init__(self, protein_data, frame_list, atom_selection="name CA"):

        protein_data.trajectory_data.trajectory[0]

        self.property = []
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
        super(frame_list)

    def sample(self, size):
        self.sampled_frame_list = random.sample(self.frame_list, size)
        super.sample(size)


class Distance:
    def __init__(self, Prop_vector, clean=False):

        min_value = np.min(Prop_vector)
        max_value = np.max(Prop_vector)

        # discretisation of the original vector with all values
        freq_prop_vector = discretize_to_dict(Prop_vector, min_value, max_value)

        if clean:
            freq_prop_vector = replace_zero(freq_prop_vector)
        # replace 0 with small number and then rescale values to make the sum
        # equal to 1
        # freq_prop_vector_clean = replace_zero(freq_prop_vector)

        # test - the distance should be zero (or close to zero) on itself
        sample_size = len(Prop_vector)

        distance = self.calculate_distance(freq_prop_vector, freq_prop_vector)

        print("size: {0:6d} distance: {1:.3f}".format(sample_size, distance))

        print("-------------------------")

        for sample_size in [10, 50, 100, 200, 400, 500]:
            # for sample_size in [800,4000,8000,16000,32000,40000]:

            sub_prop_vector = random.sample(Prop_vector, sample_size)

            # discretisation of the subsampled vector
            freq_sub_prop_vector = discretize_to_dict(
                sub_prop_vector, min_value, max_value
            )

            if clean:
                freq_sub_prop_vector = replace_zero(freq_sub_prop_vector)
            # calculate the kl divergence
            pq_distance = self.calculate_distance(
                freq_sub_prop_vector, freq_sub_prop_vector
            )

            print("size: {0:6d} distance: {1:.3f}".format(sample_size, pq_distance))

    def calculate_distance(self, x, y):
        raise NotImplementedError()


class BhattaDistance(Distance):
    def calculate_distance(self, x, y):
        return dictances.bhattacharyya(x, y)


class KLDiverDistance(Distance):
    def __init__(self, Prop_vector):
        super(Prop_vector, clean=True)

    def calculate_distance(self, x, y):
        return dictances.kullback_leibler(x, y)


class PearsonDictDistance(Distance):
    def calculate_distance(self, x, y):
        return dictances.pearson(x, y)
