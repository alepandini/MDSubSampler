import mdss_property
import mdss_dissimilarity


class DistributionDissimilarity:
    """
    A class used to compare the distributions between 2 property calculation over a trajectory

    Attributes
    ----------
    property_full :

    property_sample :

    """

    def __init__(
        self,
        property_full,
        property_sample,
        dissimilarity,
    ):
        self.property_full = property_full
        self.property_sample = property_sample
        self.dissimilarity = dissimilarity


class DistributionDissimilaritySimple(DistributionDissimilarity):
    def simple_dissimilarity_between_distributions(self):

        prop_name = self.property_full.display_name
        self.dissimilarity.dissimilarity
