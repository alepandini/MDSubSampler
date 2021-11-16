import MDSS_Protein_Sampler

# Create a ProteinData object
prot_data = MDSS_Protein_Sampler.ProteinData(
    "./Testing/MD01_1lym_example_fit_short.xtc",
    "./Testing/MD01_1lym_example.gro",
    config_parameters=None,
)

# Create a ProteinProperty object
prot_property = MDSS_Protein_Sampler.ProteinProperty(
    prot_data,
)
