import mdss_property


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
        property_full_traj,
        property_sam_traj,
        # property_class,
        # protein_data,
        # distance_class,
    ):

        self.property_full_traj = property_full_traj
        self.property_sam_traj = property_sam_traj
        self.property_full_traj = mdss_property.ProteinProperty(self.property_full_traj)
        self.property_sam_traj = mdss_property.ProteinProperty(self.property_sam_traj)
