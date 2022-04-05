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



# def compare_full_and_sample_protein(
#         property_class,
#         protein_data,
#         frame_list,
#         distance_class,
#         sampler,
#         sample_size,
#         number_of_iterations=None,
#     ):
#         print(f"Running {property_class.display_name}")
#         sampled_frame_list = sampler.sample(size)
#         # sampled_frame_list = sampler.sample(low, high, size, dtype=int)
#         prop = property_class(protein_data, frame_list, atom_selection)
#         prop.calculate_property()
#         prop_sample = property_class(protein_data, sampled_frame_list, atom_selection)
#         prop_sample.calculate_property()
#         distance_obj = distance_class(prop, prop_sample)
#         prop.write_property_vector(
#             "{}_{}_{}.dat".format(
#                 file_prefix, property_class.display_name, distance_class.display_name
#             )
#         )
#         prop_sample.write_property_vector(
#             "{}_{}_sample_{}.dat".format(
#                 file_prefix, property_class.display_name, distance_class.display_name
#             )
#         )
#         prop.write_property_discretised_vector(
#             "{}_{}_{}_{}.dat".format(
#                 file_prefix,
#                 property_class.display_name,
#                 "discr",
#                 distance_class.display_name,
#             )
#         )
#         prop_sample.write_property_discretised_vector(
#             "{}_{}_{}_sample_{}.dat".format(
#                 file_prefix,
#                 property_class.display_name,
#                 "discr",
#                 distance_class.display_name,
#             )
#         )
#         return distance_obj.distance