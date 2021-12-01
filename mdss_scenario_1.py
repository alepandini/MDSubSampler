import mdss_protein_sampler

# Need to add arg pars (python) --trajectory "name of file" --topology "name of file"
# To run it from the terminal by giving both names of files
# Get the protein data
p_data = mdss_protein_sampler.ProteinData(
    "./data/MD01_1lym_example_fit_short.xtc",
    "./data/MD01_1lym_example.gro",
    config_parameters=None,
)

# Use Random Sampling to get a sample of the protein trajectory
frame_list = list(range(1000))
prot_sample = mdss_protein_sampler.RandomSampler(frame_list, seed_number=1999)
sample = prot_sample.sample(100)

# Properties  prop
# Create a function with arguments: mdss_protein_sampler class, p_data,
# frame_list, sampled_frame_list that returns a tuple (rmsd_property, rmsd_s_property)
rmsd_property = mdss_protein_sampler.RMSDProperty(p_data, frame_list)
rmsd_s_property = mdss_protein_sampler.RMSDProperty(
    p_data, prot_sample.sampled_frame_list
)

rog_property = mdss_protein_sampler.RadiusOfGyrationProperty(p_data, frame_list)
rog_s_property = mdss_protein_sampler.RadiusOfGyrationProperty(
    p_data, prot_sample.sampled_frame_list
)

# Create function similar to property
# Calculate distance between original and sample protein
distance_rmsd = mdss_protein_sampler.Distance(rmsd_property, rmsd_s_property)
print(distance_rmsd.distance)

distance_rog = mdss_protein_sampler.Distance(rog_property, rog_s_property)
print(distance_rog.distance)
