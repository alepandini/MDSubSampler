import mdss.protein_data as pd
import parser as p
import os
import sys


def sampling_workflow(arg_list):
    """
    Implements sampling of data and measures dissimilarity between distributions

    Attributes
    -----------
    arg_list: List of arguments that were inputed by the user in the parser
    """
    args = p.parse_args(arg_list)
    print(args)

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    """
    Create property class object with user's input property
    """
    property_class = p.PROPERTY_CLASS_MAPPING[args.property]
    """
    Check if user input includes an xvg file with precalculated property
    οtherwise read protein trajectory from input files and calculate property
    """
    if args.xvg_file is not None:
        property = property_class.from_xvg(args.xvg_file)
    else:
        p_data = pd.ProteinData(
            args.trajectory_file, args.topology_file, config_parameters=None
        )
        property = property_class(p_data, args.atom_selection)
        property.calculate_property()
    """
    Create sampler class object with user's input sampler
    Depending on the sampler ensure necessary arguments were provided
    """
    sampler_class = p.SAMPLER_CLASS_MAPPING[args.sampler]
    if args.sampler == "RandomSampler":
        sampler = sampler_class(property, args.seed_number)
    if args.sampler == "UniformSampler":
        sampler = sampler_class(property, args.strata_number)
    elif args.sampler == "WeightedSampler":
        sampler = sampler_class(property, args.seed_number)
    elif args.sampler == "StratifiedSampler":
        sampler = sampler_class(property, args.layers)
    elif args.sampler == "BootstrappingSampler":
        sampler = sampler_class(property, args.number_of_iterations)
    """
    Get the sample of the vector with calculated property with user's input size
    """
    property_sample = sampler.sample(args.size)
    """
    Create dissimilarity class object with user's input dissimilarity
    """
    dissimilarity_class = p.DISSIMILARITY_CLASS_MAPPING[args.dissimilarity]
    dissimilarity_object = dissimilarity_class(property, property_sample)
    """
    Calculate dissimilarity between full and sample trajectory
    """
    dissimilarity_score = dissimilarity_object.calculate_dissimilarity()
    """
    Create file with calculated property for full trajectory
    """
    filename = "{}_{}_{}.dat".format(
        args.file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath = os.path.join(args.output_folder, filename)
    property.write_property_vector(filepath)
    """
    Create file with calculated property for sample trajectory
    """
    filename_sample = "{}_{}_sample_{}.dat".format(
        args.file_prefix, property_class.display_name, dissimilarity_class.display_name
    )
    filepath_sample = os.path.join(args.output_folder, filename_sample)
    property_sample.write_property_vector(filepath_sample)
    """
    Crete file with data report that includes important statistics about trajectory
    """
    p_data.property_data_report()
    print("Dissimilarity: {}".format(dissimilarity_score))
    return dissimilarity_score


def main(arg_list):
    sampling_workflow(arg_list)


if __name__ == "__main__":
    arg_list = sys.argv[1:]
    main(arg_list)
