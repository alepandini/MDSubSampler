import mdss_protein_sampler
import argparse

# Function that creates the parser
def parse_args():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument(
        "--traj",
        dest="trajectory_file",
        required=True,
        help="the path to trajectory file",
    )
    my_parser.add_argument(
        "--top", dest="topology_file", required=True, help="the path to topology file"
    )
    my_parser.add_argument(
        "--prefix",
        dest="file_prefix",
        required=True,
        help="the prefix for output files",
    )

    return my_parser.parse_args()


args = parse_args()
file_prefix = args.file_prefix

# Get the protein data
p_data = mdss_protein_sampler.ProteinData(
    args.trajectory_file, args.topology_file, config_parameters=None
)

# Use Bootstrapping Sampling to get a sample of the protein trajectory
frame_list = list(range(1000))
size = 100
number_of_iterations = 100
prot_sampler = mdss_protein_sampler.BootstrappingSampler(frame_list)

# Function that calculates property, distance and save in a file
def compare_full_and_sample_protein(
    property_class,
    protein_data,
    frame_list,
    distance_class,
    sampler,
    sample_size,
    number_of_iterations,
):
    print(f"Running {property_class.display_name}")
    sampled_frame_list = sampler.sample(sample_size, number_of_iterations)
    prop = property_class(protein_data, frame_list)
    prop_sample = property_class(protein_data, sampled_frame_list)
    distance_obj = distance_class(prop, prop_sample)
    prop.write_property_vector(
        "{}_{}_{}.dat".format(
            file_prefix, property_class.display_name, distance_class.display_name
        )
    )
    prop_sample.write_property_vector(
        "{}_{}_sample_{}.dat".format(
            file_prefix, property_class.display_name, distance_class.display_name
        )
    )
    # prop.write_property_discretised_vector(
    #     "{}_{}_{}_{}.dat".format(
    #         file_prefix,
    #         property_class.display_name,
    #         "discr",
    #         distance_class.display_name,
    #     )
    # )
    # prop_sample.write_property_discretised_vector(
    #     "{}_{}_{}_sample_{}.dat".format(
    #         file_prefix,
    #         property_class.display_name,
    #         "discr",
    #         distance_class.display_name,
    #     )
    # )
    return distance_obj.distance


# RMSD property and simple distance
distance_rmsd = compare_full_and_sample_protein(
    mdss_protein_sampler.RMSDProperty,
    p_data,
    frame_list,
    mdss_protein_sampler.Distance,
    prot_sampler,
    size,
    number_of_iterations,
)