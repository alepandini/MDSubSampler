import MDSS_Protein_Sampler

# Create a ProteinData object
p_data = MDSS_Protein_Sampler.ProteinData(
    "./Testing/MD01_1lym_example_fit_short.xtc",
    "./Testing/MD01_1lym_example.gro",
    config_parameters=None,
)

prot_data.output_trj_summary()
# prot_data.statistical_analysis()   Discussion with Ale about the arguments
# prot_data.add_property()           Same and about what is this method about

# Sampling the protein with the Random sampler
frame_list = list(range(1000))
prot_sample = MDSS_Protein_Sampler.RandomSampler(frame_list, seed_number=1999)
frame_list_sample = prot_sample.sample(100)

# Create a RMSDProperty object for both full protein and sample
p_property = MDSS_Protein_Sampler.RMSDProperty(p_data, frame_list)
s_property = MDSS_Protein_Sampler.RMSDProperty(p_data, frame_list_sample)

# Calculate the distance between the two proteins in terms of the RMSD property
distance = MDSS_Protein_Sampler.BhattaDistance(p_property, s_property)
print(distance.distance)
