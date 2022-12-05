import dictances
from scipy.special import rel_entr
from mdss.log_setup import log
from math import isinf


class Dissimilarity:
    """
    Represents dissimilarity measure between property of full and sample trajectory

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    n_bins : int
        number of bins for generating the discretized vector
    """

    display_name = None

    def __init__(self, target_property, ref_property, n_bins=100):
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
        Calculates dissimilarity between average values of target (sample trajectory)
        and reference (original trajectory) properties
        """
        self.dissimilarity = abs(
            self.target_property.avg_value - self.ref_property.avg_value
        )
        log.info(
            "{:15s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )
        return self.dissimilarity


class Bhattacharya(Dissimilarity):
    """
    Represents the Bhattacharya dissimilarity between target (sample traj) and reference
    (original traj) properties

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    """

    display_name = "Bhattacharyya"

    def __init__(self, target_property, ref_property):
        super().__init__(target_property, ref_property)
        self.dissimilarity_name = "Bhattacharya"
        self.dissimilarity_threshold = 0.2

    def calculate_dissimilarity(self):
        """
        Calculates Bhattacharya distance between target (sample traj) and reference
        (original traj) properties
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
    Represents the Kullback-Leibler divergence between target (sample traj) and reference
    (original traj) properties

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    """

    display_name = "KullbackLeibler"

    def __init__(self, target_property, ref_property):
        super().__init__(target_property, ref_property)
        self.dissimilarity_name = "Kullback-Leibler"
        self.dissimilarity_threshold = 0.2

    def calculate_dissimilarity(self):
        """
        Calculates Kullback-Leibler divergence between target (sample traj) and reference
        (original traj) properties
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
    Represents the Pearson distance between target (sample traj) and reference
    (original traj) properties

    Attributes
    ----------
    target_property : ProteinProperty object
        Refers to the calculated property of the sampled or target protein trajectory
    ref_property : ProteinProperty object
        Refers to the reference calculated property (e.g. for the full protein trajectory)
    """

    display_name = "Pearson"

    def __init__(self, target_property, ref_property):
        super().__init__(target_property, ref_property)
        self.dissimilarity_name = "Pearson"
        self.dissimilarity_threshold = 0.2

    def calculate_dissimilarity(self):
        """
        Calculates pearson distance between target (sample traj) and reference
        (original traj) properties
        """
        self.dissimilarity = dictances.pearson(
            self.target_property.property_distribution_dict,
            self.ref_property.property_distribution_dict,
        )
        log.info(
            "{:15s} Dissimilarity score: {:4.5f}".format("OUTPUT", self.dissimilarity)
        )
        return self.dissimilarity
