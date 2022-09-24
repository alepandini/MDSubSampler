from mdss_graph import plot_distribution
import mdss_protein_data
import mdss_property
import mdss_parser as p
from mdss_logging import logging as log
import os
import sys


def sampling_workflow(arg_list):
    args = p.parse_args(arg_list)
    print(args)

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    property_class = p.PROPERTY_CLASS_MAPPING[args.property]
    if args.xvg_file is not None:
        property = property_class.from_xvg(args.xvg_file)
    else:
        p_data = mdss_protein_data.ProteinData(
            args.trajectory_file, args.topology_file, config_parameters=None
        )
        property = property_class(p_data, args.atom_selection)
        property.calculate_property()

    sampler_class = p.SAMPLER_CLASS_MAPPING[args.sampler]
    if args.sampler == "RandomSampler":
        sampler = sampler_class(property, args.seed_number)
    if args.sampler == "UniformSampler":
        sampler = sampler_class(property, args.low, args.high, args.dtype)
    elif args.sampler == "StratifiedSampler":
        sampler = sampler_class(property, args.layers)
    elif args.sampler == "BootstrappingSampler":
        sampler = sampler_class(property, args.number_of_iterations)

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

    print(
        {
            k: {
                "min": v.min_value,
                "max": v.max_value,
                "atom_selection": v.atom_selection,
                "name": v.display_name,
            }
            for k, v in p_data.property_dict.items()
        }
    )
    return dissimilarity_score


def main(arg_list):
    sampling_workflow(arg_list)


if __name__ == "__main__":
    arg_list = sys.argv[1:]
    main(arg_list)
