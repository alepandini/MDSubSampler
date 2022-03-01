import mdss_protein_data
import mdss_property
import mdss_sampler
import mdss_distance
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
p_data = mdss_protein_data.ProteinData(
    args.trajectory_file, args.topology_file, config_parameters=None
)
# selection = p_data.frame_selection([1, 3, slice(10, 20)])
# print(selection.frames)

# Use Random Sampling to get a sample of the protein trajectory
frame_list = list(range(1000))
size = 100
prot_sampler = mdss_sampler.RandomSampler(frame_list, seed_number=1999)

# # Use Uniform Sampling to get a sample of the protein trajectory
# frame_list = list(range(1000))
# size = 100
# low = 0
# high = 1000
# prot_sampler = mdss_sampler.UniformSampler(frame_list)

# # Use Stratified Sampling to get a sample of the protein trajectory
# frame_list = list(range(1000))
# layers = [range(1, 20), range(20, 45), range(45, 100), range(100, 110)]
# size = 100
# prot_sampler = mdss_sampler.StratifiedSampler(frame_list, layers)

# # Use Bootstrapping Sampling to get a sample of the protein trajectory
# frame_list = list(range(1000))
# size = 100
# number_of_iterations = 100
# prot_sampler = mdss_sampler.BootstrappingSampler(frame_list)

# Function that calculates property, distance and save in a file
def compare_full_and_sample_protein(
    property_class,
    protein_data,
    frame_list,
    distance_class,
    sampler,
    sample_size,
    number_of_iterations=None,
):
    print(f"Running {property_class.display_name}")
    sampled_frame_list = sampler.sample(size, number_of_iterations=None)
    # sampled_frame_list = sampler.sample(low, high, size, dtype=int)
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
    prop.write_property_discretised_vector(
        "{}_{}_{}_{}.dat".format(
            file_prefix,
            property_class.display_name,
            "discr",
            distance_class.display_name,
        )
    )
    prop_sample.write_property_discretised_vector(
        "{}_{}_{}_sample_{}.dat".format(
            file_prefix,
            property_class.display_name,
            "discr",
            distance_class.display_name,
        )
    )
    return distance_obj.distance


# # RMSD property and simple distance
# distance_rmsd = compare_full_and_sample_protein(
#     mdss_property.RMSDProperty,
#     p_data,
#     list(range(1000)),
#     mdss_distance.Distance,
#     prot_sampler,
#     size,
#     number_of_iterations=None,
# )

# PCA property
pca = compare_full_and_sample_protein(
    mdss_property.PCA,
    p_data,
    list(range(1000)),
    mdss_distance.Distance,
    prot_sampler,
    size,
    number_of_iterations=None,
)

# # RMSD property and simple distance
# distance_between_atoms = compare_full_and_sample_protein(
#     mdss_property.DistanceProperty,
#     p_data,
#     list(range(1000)),
#     mdss_distance.Distance,
#     prot_sampler,
#     size,
#     number_of_iterations=None,
# )

# # RMSF property and simple distance
# distance_rmsf = compare_full_and_sample_protein(
#     mdss_protein_sampler.RMSFProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.Distance,
#     prot_sampler,
#     size,
# )

# # Radius of Gyration property and simple distance
# distance_rog = compare_full_and_sample_protein(
#     mdss_protein_sampler.RadiusOfGyrationProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.Distance,
#     prot_sampler,
#     size,
# )

# # Dihedral Angles property and simple distance
# distance_dih_ang = compare_full_and_sample_protein(
#     mdss_protein_sampler.DihedralAngles,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.Distance,
#     prot_sampler,
#     size,
# )

# # RMSD property and Bhatta distance distance
# distance_rmsd = compare_full_and_sample_protein(
#     mdss_protein_sampler.RMSDProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.BhattaDistance,
#     prot_sampler,
#     size,
# )

# # RMSF property and Bhatta distance distance
# distance_rmsf = compare_full_and_sample_protein(
#     mdss_protein_sampler.RMSFProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.BhattaDistance,
#     prot_sampler,
#     size,
# )

# # Radius of Gyration property and Bhatta distance
# distance_rog = compare_full_and_sample_protein(
#     mdss_protein_sampler.RadiusOfGyrationProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.BhattaDistance,
#     prot_sampler,
#     size,
# )

# # RMSD property and KL distance
# distance_rmsd = compare_full_and_sample_protein(
#     mdss_protein_sampler.RMSDProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.KLDiverDistance,
#     prot_sampler,
#     size,
# )

# # RMSF property and KL distance
# distance_rmsf = compare_full_and_sample_protein(
#     mdss_protein_sampler.RMSFProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.KLDiverDistance,
#     prot_sampler,
#     size,
# )


# # Radius of Gyration property and KL distance
# distance_rog = compare_full_and_sample_protein(
#     mdss_protein_sampler.RadiusOfGyrationProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.KLDiverDistance,
#     prot_sampler,
#     size,
# )
# # RMSD property and Pearson distance
# distance_rmsd = compare_full_and_sample_protein(
#     mdss_protein_sampler.RMSDProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.PearsonDictDistance,
#     prot_sampler,
#     size,
# )

# # RMSF property and Pearson distance
# distance_rmsf = compare_full_and_sample_protein(
#     mdss_protein_sampler.RMSFProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.PearsonDictDistance,
#     prot_sampler,
#     size,
# )

# # Radius of Gyration property and Pearson distance
# distance_rog = compare_full_and_sample_protein(
#     mdss_protein_sampler.RadiusOfGyrationProperty,
#     p_data,
#     list(range(1000)),
#     mdss_protein_sampler.PearsonDictDistance,
#     prot_sampler,
#     size,
# )
