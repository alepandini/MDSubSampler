import dictances
from mdss_logging import logging as log


class Dissimilarity:
    """
    A class used to calculate the dissimilarity distance in terms of property between
    a full protein trajectory and a sample of it to identify the difference

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    """

    display_name = None

    def __init__(self, target_property, ref_property):
        self.target_property = target_property
        self.ref_property = ref_property
        self.min_value = min(
            min(target_property.property_vector), min(ref_property.property_vector)
        )
        self.max_value = max(
            max(target_property.property_vector), max(ref_property.property_vector)
        )
        self.target_property.discretize_vector(min_value=self.min_value, max_value=self.max_value)
        self.ref_property.discretize_vector(min_value=self.min_value, max_value=self.max_value)
        self.dissimilarity = None

    def calculate_dissimilarity(self):
        """
        Method that calculates the difference between the average values of the
        two calculated property vectors.
        """
        self.dissimilarity = self.target_property.avg_value - self.ref_property.avg_value
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )


class BhattaCoefficient(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Bhattacharya coefficient between
    two property vectors

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    """

    display_name = "bhatta"

    def __init__(self, target_property, ref_property):
        super().__init__(target_property, ref_property)

    def calculate_dissimilarity(self):
        """
        Method that returns the Bhatta coefficient from distance between two vectors
        """
        self.dissimilarity = dictances.bhattacharyya(
            self.target_property.property_distribution_dict,
            self.ref_property.property_distribution_dict
        )
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )


class KLDivergence(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Kullback-Leibler divergence between
    two property vectors

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    """

    display_name = "kl"

    def __init__(self, target_property, ref_property):
        super().__init__(target_property, ref_property)

    def calculate_dissimilarity(self):
        """
        Method that returns the KL divergence distance between two vectors
        """
        self.dissimilarity = dictances.kullback_leibler(
            self.target_property.property_distribution_dict,
            self.ref_property.property_distribution_dict,
        )
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )


class PearsonCoefficient(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Pearson coefficient between
    two property vectors

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    """

    display_name = "pearson"

    def __init__(self, target_property, ref_property):
        super().__init__(target_property, ref_property)

    def calculate_dissimilarity(self):
        """
        Method that returns the pearson coefficient between two vectors
        """
        self.dissimilarity = dictances.pearson(
            self.target_property.property_distribution_dict,
            self.ref_property.property_distribution_dict,
        )
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )
