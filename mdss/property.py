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
import numpy as np
import MDAnalysis as mda
from mdss.protein_data import ProteinData
from mdss.log_setup import log
from mdss.dissimilarity import *


class ProteinProperty:
    """
    Class representing protein property on a single or multiple frames of a protein trajectory.

    Attributes
    ----------
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    atom_selection : str, optional
        Atom selection for property calculation. Default string value is CA atoms"
    """

    display_name = None

    def __init__(self, protein_data, atom_selection="name CA"):
        if not isinstance(protein_data, ProteinData):
            log.warning(
                "{:12s} An instance of ProteinData is required. protein_data attribute set to None".format(
                    "STEPS"
                )
            )
            protein_data = None
        self.protein_data = protein_data
        self.atom_selection = atom_selection
        self.ref_coordinates = []
        self.ref_frame_index = None
        self.property_vector = []
        self.discretized_property_vector = []
        self.property_distribution_dict = {}
        self.frame_indices = []
        self.ref_dissimilarity = 0
        self.recalculate = False
        self.property_key = self._add_reference_to_protein_data()
        log.info("{:15s} Atom selection: {}".format("INPUT", self.atom_selection))
        log.info("{:15s} Property name: {}".format("INPUT", self.display_name))

    @classmethod
    def from_xvg(cls, xvg_filepath):
        """
        Class method that constructs an instance of the class directly from an xvg file.

        Parameters
        ----------

        xvg_filepath : str
            Path where the xvg file is located.

        Returns
        -------
        object
            An instance of the class that was used to call the method.
        """
        instance = cls(protein_data=None, atom_selection=[None, None])
        _frames, distance_values = np.loadtxt(xvg_filepath, unpack=True)
        instance.property_vector = distance_values
        instance._property_statistics()
        instance.discretize_vector()
        return instance

    def _add_reference_to_protein_data(self):
        """
        Link ProteinProperty to ProteinData class.

        Returns
        -------
        string
            Display name of reference property with unique identifier.

        """
        if self.protein_data is not None:
            property_key = self.protein_data.add_property(self, self.display_name)
            return property_key

    def discretize_vector(self, min_value=None, max_value=None, n_bins=100):
        """
        Discretise vector for Bhattacharyya and Kullback-Leibler distances.

        Parameters
        ----------
        min_value : int
            Minimum value for vector.
        max_value : int
            Maximum value for vector.
        n_bins : int, optional
            Number of bins for discretized vector. Default value is 100.
        """
        if min_value is None:
            min_value = self.min_value
        if max_value is None:
            max_value = self.max_value
        bin_size = (max_value - min_value) / n_bins
        bin_vector = np.arange(min_value, max_value, bin_size)
        self.discretized_property_vector = [
            bin_vector[i - 1] for i in np.digitize(self.property_vector, bin_vector)
        ]
        counts, bins = np.histogram(self.property_vector, bins=bin_vector)
        self.property_distribution_dict = dict(
            zip(bins, (counts / len(self.property_vector)))
        )

    def _property_statistics(self):
        """
        Calculate minimum, maximum and average values for a specific vector.
        """
        self.min_value = np.min(self.property_vector)
        self.max_value = np.max(self.property_vector)
        self.avg_value = np.average(self.property_vector)

    def set_reference_coordinates(self, frame_index=None):
        """
        Set up on a reference frame and extract copy of coordinates.

        Parameters
        ----------
        frame_index : int
            Reference structure (i.e. frame) from inputed protein trajectory.

        Returns
        -------
        Boolean
            False if property is not associated to a protein data object, otherwise True.
        """
        if self.protein_data is not None:
            if frame_index is None:
                trj_data = self.protein_data.topology_data
            else:
                trj_data = self.protein_data.trajectory_data
                trj_data.trajectory[frame_index]
            if isinstance(self.atom_selection, list):
                self.ref_coordinates = [
                    trj_data.select_atoms(selection).positions.copy()
                    for selection in self.atom_selection
                ]
            else:
                self.ref_coordinates = trj_data.select_atoms(
                    self.atom_selection
                ).positions.copy()
                self.ref_frame_index = frame_index
            return True
        else:
            print("Warning: property is not associated to a protein data object.")
            log.warning(
                "{:12s} Property is not associated to a protein data object".format(
                    "STEPS"
                )
            )
            return False

    def write_property_vector(self, outfilepath):
        """
        Generate and save file with calculation of reference property for all frames in protein trajectory.

        Parameters
        ----------
        outfilepath : str
            Path where output file is saved.
        """
        if outfilepath is not None:
            with open(outfilepath, "w") as f:
                for i, value in zip(self.frame_indices, self.property_vector):
                    f.write("{} {}\n".format(i, value))

    def write_discretized_property_vector(self, outfilepath):
        """
        Generate and save file with calculation of discretised reference property for all frames in protein trajectory.

        Parameters
        ----------
        outfilepath : str
            Path where output file is saved.
        """
        if len(self.property_distribution_dict) < 1:
            self.discretize_vector()
        if outfilepath is not None:
            with open(outfilepath, "w") as f:
                for i, value in zip(
                    self.frame_indices, self.discretized_property_vector
                ):
                    f.write("{} {}\n".format(i, value))

    def write_property_distribution_dict(self, outfilepath):
        """
        Generate and save dictionary file which contains information for property distribution values.

        Parameters
        ----------
        outfilepath : str
            Path where output file is saved.
        """
        if len(self.property_distribution_dict) < 1:
            self.discretize_vector()
        if outfilepath is not None:
            with open(outfilepath, "w") as f:
                for key, value in self.property_distribution_dict.items():
                    f.write("{} {}\n".format(key, value))


class SampledProperty(ProteinProperty):
    """
    Subclass of ProteinProperty class representing sampled protein property on individual frames from trajectory.

    Attributes
    ----------
    original_property : vector
        Contains calculated reference property for full trajectory.
    sampled_property_vector : vector
        Contains calculated reference property for sampled trajectory.
    sampled_frame_indices : list
        Contains frame indices for sampled trajectory.
    samples_indices : list
        Contains samples indices.
    sample_size : int
        Size of sampled trajectory.
    dissimilarity_measure : Dissimilarity
        An instance of the Dissimilarity class representing the dissimilarity measure.
    """

    def __init__(
        self,
        original_property,
        sampled_property_vector,
        sampled_frame_indices,
        samples_indices,
        sample_size,
        dissimilarity_measure=Bhattacharyya,
    ):
        self.protein_data = original_property.protein_data
        self.atom_selection = original_property.atom_selection
        self.ref_property = original_property
        self.ref_coordinates = []
        self.property_vector = sampled_property_vector
        self.discretized_property_vector = []
        self.property_distribution_dict = {}
        self.frame_indices = sampled_frame_indices
        self.samples_indices = samples_indices
        self.dissimilarity_name = None
        self.dissimilarity_threshold = None
        self._property_statistics()
        self.ref_dissimilarity = self._get_dissimilarity_to_ref(dissimilarity_measure)
        self.display_name = "{}_{}".format(sample_size, original_property.display_name)
        self.property_key = self._add_reference_to_protein_data()

    def _get_dissimilarity_to_ref(self, dissimilarity_measure):
        """
        Measure dissimilarity between full and sampled calculatd reference property.

        Parameters
        ----------
        dissimilarity_measure : Dissimilarity
            An instance of the Dissimilarity class representing the dissimilarity measure.

        Returns
        -------
        Dissimilarity
            An instance of Dissimilarity class.
        """
        diss_obj = dissimilarity_measure(self, self.ref_property)
        diss_obj.calculate_dissimilarity()
        self.dissimilarity_name = diss_obj.dissimilarity_name
        self.dissimilarity_threshold = diss_obj.dissimilarity_threshold
        return diss_obj.dissimilarity

    def calculate_property(self):
        pass

    def get_samples_averages(self):
        """
        Retrieve average value for each sample trajectory.

        Returns
        -------
        dictionary
            Dictionary with sample sizes and their average values.

        """
        samples_labels = set(self.samples_indices)
        average_dict = {}
        for label in samples_labels:
            average_value = np.mean(
                [
                    value
                    for value, i in zip(self.property_vector, self.samples_indices)
                    if i == label
                ]
            )
            average_dict[label] = average_value
        return average_dict

    def get_average(self):
        """
        Retrieve average value from samples averages.

        Returns
        -------
        float
            Average value among all values in a dictionary.
        """
        average_dict = self.get_samples_averages()
        average_value = np.mean(list(average_dict.values()))
        return average_value
