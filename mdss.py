from mdss_graph import plot_distribution
import mdss_protein_data
import mdss_property
import mdss_parser as p
from mdss_logging import logging as log
import os
import sys


# def run_subsampler(p_data, property_class, sampler_class):
#     """
#     Method that uses the user input for property calculation, sampling method
#     and size of sample and returns a subsample trajectory along with a log
#     file with diagnostics for the particular property."""
#     log.logging.info("The MDSubsampler is running:")
#     property.calculate_property()
#     property_sample.calculate_property()
#     print(f"Calculating {property_class.display_name}")
#     print(f"Applying {property_class.display_name}")
#     property.calculate_property()


def sampling_workflow(arg_list):
    args = p.parse_args(arg_list)
    print(args)

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    p_data = mdss_protein_data.ProteinData(
        args.trajectory_file, args.topology_file, config_parameters=None
    )

    property_class = p.PROPERTY_CLASS_MAPPING[args.property]
    # property = property_class(p_data, args.atom_selection)
    property = property_class.from_xvg("./data/ApoADK_protein_05_G55-P127.xvg")
    # property.calculate_property()

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
    sampled_indices = sampler.sampled_indices
    property_sample = mdss_property.SampledProperty(
        sampled_property_vector, sampled_indices
    )

    dissimilarity_class = p.DISSIMILARITY_CLASS_MAPPING[args.dissimilarity]
    dissimilarity = dissimilarity_class(property, property_sample)
    dissimilarity_score = dissimilarity.calculate_dissimilarity()
    print("Dissimilarity: {}".format(dissimilarity_score))
    log.info("Dissimilarity: {}".format(dissimilarity_score))

    filename = "{}_{}_{}.dat".format(
        args.file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath = os.path.join(args.output_folder, filename)

    property.write_property_vector(filepath)

    filename_sample = "{}_{}_sample_{}.dat".format(
        args.file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath_sample = os.path.join(args.output_folder, filename_sample)

    property_sample.write_property_vector_sample(filepath_sample)

    plot_distribution(
        property,
        filepath,
        "{}_{}_{}".format(
            args.file_prefix,
            property_class.display_name,
            dissimilarity_class.display_name,
        ),
        args.output_folder,
    )

    plot_distribution(
        property_sample,
        filepath_sample,
        "{}_{}_sample_{}".format(
            args.file_prefix,
            property_class.display_name,
            dissimilarity_class.display_name,
        ),
        args.output_folder,
    )

    return dissimilarity_score


def main(arg_list):
    sampling_workflow(arg_list)


if __name__ == "__main__":
    arg_list = sys.argv[1:]
    main(arg_list)

# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "11%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='BootstrappingSampler' --n-iterations=50 --size=11000 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "0.1%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1991 --size=100 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "0.5%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1992 --size=500 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "1%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1993 --size=1000 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "5%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1994 --size=5000 --dissimilarity='BhattaCoefficient'
# python mdss.py --traj "data/user.xtc" --top "data/user.gro" --prefix "10%" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='RandomSampler' --seed-number=1995 --size=10000 --dissimilarity='BhattaCoefficient'

# 0.05: 0.4847413787731183
# 0.1:0.1892646854226224
# 0.5:0.03408998524296803
# 1: 0.008814553862456592
# 5:0.003142312466437116
# 10:0.001085055552200252
