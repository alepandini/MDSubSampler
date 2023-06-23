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
import dictances
from scipy.special import rel_entr
from mdss.log_setup import log
from math import isinf


class Dissimilarity:
    """
    Class representing dissimilarity measure of calculated property between full and sample trajectory.

    Attributes
    ----------
    target_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference property for sampled trajectory.
    ref_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference property for full trajectory.
    n_bins : int, optional
        Number of bins for generation the discretized vector. Default value is 100.
    """

    display_name = None

    def __init__(self, target_property, ref_property, n_bins=100):
        """
        Initialize the Dissimilarity object.

        Parameters
        ----------
        target_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference property for the sampled trajectory.
        ref_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference property for the full trajectory.
        n_bins : int, optional
            Number of bins for generating the discretized vector. Default value is 100.
        """
        self.dissimilarity_name = "average"
        self.target_property = target_property
        self.ref_property = ref_property
        self.min_value = min(
            min(target_property.property_vector), min(ref_property.property_vector)
        )
        self.max_value = max(
            max(target_property.property_vector), max(ref_property.property_vector)
        )
        self.target_property.discretize_vector(
            min_value=self.min_value, max_value=self.max_value, n_bins=n_bins
        )
        self.ref_property.discretize_vector(
            min_value=self.min_value, max_value=self.max_value, n_bins=n_bins
        )
        self.dissimilarity = None
        self.dissimilarity_threshold = (self.max_value - self.min_value) / 0.1

    def calculate_dissimilarity(self):
        """
        Calculate dissimilarity between average property value of sample and full trajectory.

        Returns
        -------
        float
            The difference between the two average property values.
        """
        self.dissimilarity = abs(
            self.target_property.avg_value - self.ref_property.avg_value
        )
        log.info(
            "{:15s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )
        return self.dissimilarity


class Bhattacharyya(Dissimilarity):
    """
    Subclass of Dissimilarity class that represents Bhattacharyya dissimilarity measure.
    It calculates the distance of property distribution between full and sample trajectory.

    Attributes
    ----------
    target_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference
        property for sampled trajectory.
    ref_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference
        property for full trajectory.
    """

    display_name = "Bhattacharyya"

    def __init__(self, target_property, ref_property):
        """
        Initialize the Bhattacharyya object.

        Parameters
        ----------
        target_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference property for the sampled trajectory.
        ref_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference property for the full trajectory.
        """
        super().__init__(target_property, ref_property)
        self.dissimilarity_name = "Bhattacharyya"
        self.dissimilarity_threshold = 0.2

    def calculate_dissimilarity(self):
        """
        Calculate Bhattacharyya distance between average property value of sample and full trajectory.

        Returns
        -------
        float
            The Bhattacharyya distance between the two average property values.
        """
        self.dissimilarity = dictances.bhattacharyya(
            self.target_property.property_distribution_dict,
            self.ref_property.property_distribution_dict,
        )

        log.info(
            "{:15s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )
        return self.dissimilarity


class KullbackLeibler(Dissimilarity):
    """
    Subclass of Dissimilarity class that represents Kullback-Leibler divergence.
    It calculates the distance of property distribution between full and sample trajectory.

    Attributes
    ----------
    target_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference
        property for sampled trajectory.
    ref_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference
        property for full trajectory.
    """

    display_name = "KullbackLeibler"

    def __init__(self, target_property, ref_property):
        """
        Initialize the KullbackLeibler object.

        Parameters
        ----------
        target_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference
            property for the sampled trajectory.
        ref_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference
            property for the full trajectory.
        """
        super().__init__(target_property, ref_property)
        self.dissimilarity_name = "Kullback-Leibler"
        self.dissimilarity_threshold = 0.2

    def calculate_dissimilarity(self):
        """.
        Calculate Kullback-Leibler divergence between average property value of sample and full trajectory.

        Returns
        -------
        float
            The Kullback-Leibler distance between the two average property values.
        """
        P = list(self.target_property.property_distribution_dict.values())
        Q = list(self.ref_property.property_distribution_dict.values())
        rel_entropy_vector = rel_entr(P, Q)
        self.dissimilarity = sum([v for v in rel_entropy_vector if not isinf(v)])
        log.info(
            "{:15s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )
        return self.dissimilarity


class Pearson(Dissimilarity):
    """
    Subclass of Dissimilarity class that represents Pearson distance.
    It calculates the distance of property distribution between full and sample trajectory.

    Attributes
    ----------
    target_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference
        property for sampled trajectory.
    ref_property : ProteinProperty
        An instance of the ProteinProperty class representing the calculated reference
        property for full trajectory.
    """

    display_name = "Pearson"

    def __init__(self, target_property, ref_property):
        """
        Initialize the Pearson object.

        Parameters
        ----------
        target_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference
            property for the sampled trajectory.
        ref_property : ProteinProperty
            An instance of the ProteinProperty class representing the calculated reference
            property for the full trajectory.
        """
        super().__init__(target_property, ref_property)
        self.dissimilarity_name = "Pearson"
        self.dissimilarity_threshold = 0.2

    def calculate_dissimilarity(self):
        """
        Calculate Pearson distance between average property value of sample and full trajectory.

        Returns
        -------
        float
            The Pearson distance between the two average property values.
        """
        self.dissimilarity = dictances.pearson(
            self.target_property.property_distribution_dict,
            self.ref_property.property_distribution_dict,
        )
        log.info(
            "{:15s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )
        return self.dissimilarity


"""
    Dictionary mapping dissimilarity names to their corresponding dissimilarity classes.

    Keys:
    -----
    Bhattacharyya : Bhattacharyya
        Bhattacharyya dissimilarity class.
    KullbackLeibler : KullbackLeibler
        Kullback-Leibler dissimilarity class.
    Pearson : Pearson
        Pearson dissimilarity class.
    Dissimilarity : Dissimilarity
        Generic dissimilarity class.

"""
dissimilarity_class_dict = {
    "Bhattacharyya": Bhattacharyya,
    "KullbackLeibler": KullbackLeibler,
    "Pearson": Pearson,
    "Dissimilarity": Dissimilarity,
}
