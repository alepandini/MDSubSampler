from numpy import require
import mdss_protein_data
import mdss_property
import mdss_parser as p
import operator as op
import os

if __name__ == "__main__":
    args = p.parse_args()
    print(args)

    file_prefix = args.file_prefix
    output_folder = args.output_folder
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    p_data = mdss_protein_data.ProteinData(
        args.trajectory_file, args.topology_file, config_parameters=None
    )

    property_class = p.PROPERTY_CLASS_MAPPING[args.property]
    property = property_class(p_data, args.atom_selection)
    property.calculate_property()

    sampler_class = p.SAMPLER_CLASS_MAPPING[args.sampler]
    if args.sampler == "RandomSampler":
        sampler = sampler_class(property.property_vector, args.seed_number)
    if args.sampler == "UniformSampler":
        sampler = sampler_class(
            property.property_vector, args.low, args.high, args.dtype
        )
    elif args.sampler == "StratifiedSampler":
        sampler = sampler_class(property.property_vector, args.layers)
    elif args.sampler == "BootstrappingSampler":
        sampler = sampler_class(property.property_vector, args.number_of_iterations)

    sampled_property_vector = sampler.sample(args.size)
    property_sample = mdss_property.SampledProperty(sampled_property_vector)

    dissimilarity_class = p.DISSIMILARITY_CLASS_MAPPING[args.dissimilarity]
    dissimilarity = dissimilarity_class(property, property_sample)
    print("Dissimilarity: {}".format(dissimilarity.calculate_dissimilarity()))
    distribution = p.mdss_distribution.DistributionDissimilaritySimple(
        property, property_sample, dissimilarity
    )
    distribution.simple_dissimilarity_between_distributions()

    distrib = p.mdss_distribution.DistributionDissimilarity(
        property, property_sample, dissimilarity
    )
    filename = "{}_{}_{}.dat".format(
        file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath = os.path.join(args.output_folder, filename)

    property.write_property_vector(filepath)

    filename_sample = "{}_{}_sample_{}.dat".format(
        file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath_sample = os.path.join(args.output_folder, filename_sample)

    property_sample.write_property_vector(filepath_sample)

    filename_plot = "{}_{}_{}_plot.png".format(
        file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath_plot = os.path.join(args.output_folder, filename_plot)
    property.plot_property(filepath, filepath_plot)

    filename_sample_plot = "{}_{}_sample_{}_plot.png".format(
        file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath_sample_plot = os.path.join(args.output_folder, filename_sample_plot)
    property_sample.plot_property(filepath_sample, filepath_sample_plot)


""""
    Method that uses the user input for property calculation, sampling method
    and size of sample and returns a subsample trajectory along with a log
    file with diagnostics for the particular property.
    """


def run_subsampler(p_data, property_class, sampler_class):
    property.calculate_property()
    property_sample.calculate_property()
    print(
        f"Calculating {property_class.display_name}"
    )  # Will be replaced with log file
    print(f"Applying {property_class.display_name}")
    property.calculate_property()


# python mdss.py --traj "data/MD01_1lym_example_fit_short.xtc" --top "data/MD01_1lym_example.gro" --prefix "001" --output-folder "data/results" --property='RMSDProperty' --atom-selection='name CA' --sampler='RandomSampler' --seed-number=1999 --size=100 --distance='Distance'
