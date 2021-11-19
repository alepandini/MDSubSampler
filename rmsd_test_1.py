import mdss_protein_sampler

# Create a ProteinData object
p_data = mdss_protein_sampler.ProteinData(
    "./data/MD01_1lym_example_fit_short.xtc",
    "./data/MD01_1lym_example.gro",
    config_parameters=None,
)

#prot_data.output_trj_summary()
# prot_data.statistical_analysis()   Discussion with Ale about the arguments
# prot_data.add_property()           Same and about what is this method about

# Sampling the protein with the Random sampler
frame_list = list(range(1000))
prot_sample = mdss_protein_sampler.RandomSampler(frame_list, seed_number=1999)
prot_sample.sample(100)

# Create a RMSDProperty object for both full protein and sample
p_property = mdss_protein_sampler.RMSDProperty(p_data, frame_list)
s_property = mdss_protein_sampler.RMSDProperty(p_data, prot_sample.sampled_frame_list)

# Calculate the distance between the two proteins in terms of the RMSD property
distance = mdss_protein_sampler.Distance(p_property, s_property)
print(distance.distance)
