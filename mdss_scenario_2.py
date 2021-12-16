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

# Use Stratified Sampling to get a sample of the protein trajectory
frame_list = list(range(1000))
layers = [range(1, 20), range(20, 45), range(45, 100), range(100, 110)]
size = 40
prot_sampler = mdss_protein_sampler.StratifiedSampler(frame_list, layers)

# Function that calculates property, distance and save in a file
def compare_full_and_sample_protein(
    property_class,
    prop_name,
    protein_data,
    frame_list,
    distance_class,
    sampler,
    sample_size,
):
    sampled_frame_list = sampler.sample(sample_size)
    prop = property_class(protein_data, frame_list)
    prop_sample = property_class(protein_data, sampled_frame_list)
    distance_obj = distance_class(prop, prop_sample)
    prop.write_property_vector("{}_{}.dat".format(file_prefix, prop_name))
    prop_sample.write_property_vector("{}_{}_sample.dat".format(file_prefix, prop_name))
    prop.write_property_discretised_vector(
        "{}_{}_{}.dat".format(file_prefix, prop_name, "discr")
    )
    prop_sample.write_property_discretised_vector(
        "{}_{}_{}_sample.dat".format(file_prefix, prop_name, "discr")
    )
    return distance_obj.distance


# RMSD property and simple distance
distance_rmsd = compare_full_and_sample_protein(
    mdss_protein_sampler.RMSDProperty,
    "rmsd",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.Distance,
    prot_sampler,
    size,
)

# Radius of Gyration property and simple distance
distance_rog = compare_full_and_sample_protein(
    mdss_protein_sampler.RadiusOfGyrationProperty,
    "rog",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.Distance,
    prot_sampler,
    size,
)


# RMSD property and Bhatta distance distance
distance_rmsd = compare_full_and_sample_protein(
    mdss_protein_sampler.RMSDProperty,
    "rmsd",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.BhattaDistance,
    prot_sampler,
    size,
)

# Radius of Gyration property and Bhatta distance
distance_rog = compare_full_and_sample_protein(
    mdss_protein_sampler.RadiusOfGyrationProperty,
    "rog",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.BhattaDistance,
    prot_sampler,
    size,
)

# RMSD property and KL distance
distance_rmsd = compare_full_and_sample_protein(
    mdss_protein_sampler.RMSDProperty,
    "rmsd",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.KLDiverDistance,
    prot_sampler,
    size,
)

# Radius of Gyration property and KL distance
distance_rog = compare_full_and_sample_protein(
    mdss_protein_sampler.RadiusOfGyrationProperty,
    "rog",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.KLDiverDistance,
    prot_sampler,
    size,
)
# RMSD property and Pearson distance
distance_rmsd = compare_full_and_sample_protein(
    mdss_protein_sampler.RMSDProperty,
    "rmsd",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.PearsonDictDistance,
    prot_sampler,
    size,
)

# Radius of Gyration property and Pearson distance
distance_rog = compare_full_and_sample_protein(
    mdss_protein_sampler.RadiusOfGyrationProperty,
    "rog",
    p_data,
    list(range(1000)),
    mdss_protein_sampler.PearsonDictDistance,
    prot_sampler,
    size,
)
