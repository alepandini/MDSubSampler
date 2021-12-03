import mdss_protein_sampler
import argparse

# Create and call parser
my_parser = argparse.ArgumentParser()
my_parser.add_argument(
    "--traj", dest="trajectory_file", required=True, help="the path to trajectory file"
)
my_parser.add_argument(
    "--top", dest="topology_file", required=True, help="the path to topology file"
)

args = my_parser.parse_args()

# Get the protein data
p_data = mdss_protein_sampler.ProteinData(
    args.trajectory_file, args.topology_file, config_parameters=None
)

# Use Random Sampling to get a sample of the protein trajectory
frame_list = list(range(1000))
prot_sample = mdss_protein_sampler.RandomSampler(frame_list, seed_number=1999)
sample = prot_sample.sample(100)

# Function that calculates the property of the full and the sample protein
def calculate_property(property_class, protein_data, frame_list, sampled_frame_list):
    prop = property_class(protein_data, frame_list)
    prop_sample = property_class(protein_data, sampled_frame_list)
    return (prop, prop_sample)


(rmsd_property, rmsd_sample_property) = calculate_property(
    mdss_protein_sampler.RMSDProperty,
    p_data,
    frame_list,
    prot_sample.sampled_frame_list,
)

(rog_property, rog_sample_property) = calculate_property(
    mdss_protein_sampler.RadiusOfGyrationProperty,
    p_data,
    frame_list,
    prot_sample.sampled_frame_list,
)

# Function that calculates the difference between the full trajectory and the sample
# trajectory of the protein
def calculate_distance(distance_class, property, sample_property):
    distance_obj = distance_class(property, sample_property)
    return distance_obj.distance


# Calculate simple distance for all properties
distance_rmsd = calculate_distance(
    mdss_protein_sampler.Distance, rmsd_property, rmsd_sample_property
)

distance_rog = calculate_distance(
    mdss_protein_sampler.Distance, rog_property, rog_sample_property
)

# Calculate Bhatta distance for all properties
bhatta_distance_rmsd = calculate_distance(
    mdss_protein_sampler.BhattaDistance, rmsd_property, rmsd_sample_property
)

bhatta_distance_rog = calculate_distance(
    mdss_protein_sampler.BhattaDistance, rog_property, rog_sample_property
)


# Calculate KL distance for all properties
kl_distance_rmsd = calculate_distance(
    mdss_protein_sampler.KLDiverDistance, rmsd_property, rmsd_sample_property
)

kl_distance_rog = calculate_distance(
    mdss_protein_sampler.KLDiverDistance, rog_property, rog_sample_property
)

# Calculate Pearson distance for all properties
pears_distance_rmsd = calculate_distance(
    mdss_protein_sampler.PearsonDictDistance, rmsd_property, rmsd_sample_property
)

pears_distance_rog = calculate_distance(
    mdss_protein_sampler.PearsonDictDistance, rog_property, rog_sample_property
)
print(distance_rmsd)
print(bhatta_distance_rmsd)
print(kl_distance_rmsd)
print(pears_distance_rmsd)
