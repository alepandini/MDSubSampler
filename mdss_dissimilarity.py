import dictances
from mdss_logging import logging as log


class Dissimilarity:
    """
    A class used to calculate the dissimilarity distance in terms of property between
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
        self.dissimilarity = self.calculate_dissimilarity()

    def calculate_dissimilarity(self):
        """
        Method that calculates the difference between the average values of the
        two calculated property vectors.
        """
        dissimilarity_score = self.property_1.avg_value - self.property_2.avg_value
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", dissimilarity_score)
        )
        return self.property_1.avg_value - self.property_2.avg_value


class BhattaCoefficient(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Bhattacharya coefficient between
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

    def calculate_dissimilarity(self):
        """
        Method that returns the Bhatta coefficient from distance between two vectors
        """
        property_1_discretized = self.property_1.discretize_vector(
            min_value=self.min_value, max_value=self.max_value
        )
        property_2_discretized = self.property_2.discretize_vector(
            min_value=self.min_value, max_value=self.max_value
        )
        dissimilarity_score = dictances.bhattacharyya(
            property_1_discretized, property_2_discretized
        )
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", dissimilarity_score)
        )

        return dictances.bhattacharyya(property_1_discretized, property_2_discretized)


class KLDivergence(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Kullback-Leibler divergence between
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

    def calculate_dissimilarity(self):
        """
        Method that returns the KL divergence distance between two vectors
        """
        dissimilarity_score = dictances.kullback_leibler(
            self.property_1.property_distribution_dict,
            self.property_2.property_distribution_dict,
        )
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", dissimilarity_score)
        )
        return dictances.kullback_leibler(
            self.property_1.property_distribution_dict,
            self.property_2.property_distribution_dict,
        )


class PearsonCoefficient(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Pearson coefficient between
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

    def calculate_dissimilarity(self):
        """
        Method that returns the pearson coefficient between two vectors
        """
        dissimilarity_score = dictances.pearson(
            self.property_1.property_distribution_dict,
            self.property_2.property_distribution_dict,
        )
        log.info(
            "{:18s} Dissimilarity score: {:4.5f}".format("OUTPUT", dissimilarity_score)
        )
        return dictances.pearson(
            self.property_1.property_distribution_dict,
            self.property_2.property_distribution_dict,
        )
