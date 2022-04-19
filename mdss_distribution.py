import mdss_property
import mdss_distance


class DistributionDistance:
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
        distance,
    ):
        self.property_full = property_full
        self.property_sample = property_sample
        self.distance = distance


class DistributionDistanceSimple(DistributionDistance):
    def simple_distance_between_distributions(self):

        prop_name = self.property_full.display_name
        self.distance.distance
