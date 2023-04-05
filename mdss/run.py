import mdss.protein_data as pd
import mdss.parser as p
from mdss.dissimilarity import *
from mdss.log_setup import log
import mdss.graph as g
from mdss.utilities import write_output_files
from mdss.utilities import plot_property

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

    log.info("{:15s} Prefix: {:4.5s}".format("INPUT", args.file_prefix))
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    """
    Create property class object with user's input property
    """
    property_class = p.PROPERTY_CLASS_MAPPING[args.property]
    """
    Check if user input includes an xvg file with precalculated property
    otherwise read protein trajectory from input files and calculate property
    """
    if args.xvg_file is not None:
        property = property_class.from_xvg(args.xvg_file)
    else:
        p_data = pd.ProteinData(
            args.trajectory_file,
            args.topology_file,
            config_parameters=None,
        )
        property = property_class(p_data, args.atom_selection, args.fit)
        property.calculate_property()
    """
    Create sampler class object with user's input sampler
    Depending on the sampler ensure necessary arguments were provided
    """
    sampler_class = p.SAMPLER_CLASS_MAPPING[args.sampler]
    if args.sampler == "RandomSampler":
        sampler = sampler_class(
            protein_property=property,
            protein_data=p_data,
            seed_number=args.seed_number,
            output_folder=args.output_folder,
            file_prefix=args.file_prefix,
            dissimilarity_measure=dissimilarity_class_dict[args.dissimilarity],
        )
    if args.sampler == "UniformSampler":
        sampler = sampler_class(
            protein_property=property,
            protein_data=p_data,
            strata_number=args.strata_number,
            output_folder=args.output_folder,
            file_prefix=args.file_prefix,
            dissimilarity_measure=dissimilarity_class_dict[args.dissimilarity],
        )
    elif args.sampler == "WeightedSampler":
        sampler = sampler_class(
            protein_property=property,
            protein_data=p_data,
            seed_number=args.seed_number,
            output_folder=args.output_folder,
            file_prefix=args.file_prefix,
            weights_vector=args.weights_vector,
            dissimilarity_measure=dissimilarity_class_dict[args.dissimilarity],
        )
    elif args.sampler == "StratifiedSampler":
        sampler = sampler_class(
            protein_property=property,
            protein_data=p_data,
            output_folder=args.output_folder,
            file_prefix=args.file_prefix,
            strata_vector=args.strata_vector,
            dissimilarity_measure=dissimilarity_class_dict[args.dissimilarity],
        )
    elif args.sampler == "BootstrappingSampler":
        sampler = sampler_class(
            protein_property=property,
            protein_data=p_data,
            output_folder=args.output_folder,
            file_prefix=args.file_prefix,
            number_of_iterations=args.number_of_iterations,
            seed_number=args.seed_number,
            dissimilarity_measure=dissimilarity_class_dict[args.dissimilarity],
        )
    """
    Generate all output files in case of list of sample sizes as input
    """
    if isinstance(args.size, list):
        sampled_property = sampler.scan_sample_size(
            perc_vector=args.size,
            dissimilarity_threshold=None,
            step_recording=args.step_recording,
        )
    else:
        sampled_property = sampler.sample(round(int(args.size) * p_data.n_frames / 100))

    """
    Create dissimilarity class object with user's input dissimilarity
    """
    dissimilarity_class = p.DISSIMILARITY_CLASS_MAPPING[args.dissimilarity]
    dissimilarity_object = dissimilarity_class(property, sampled_property)
    """
    Generate file with calculated property for full trajectory
    """
    filename = "{}_{}.dat".format(args.file_prefix, property_class.display_name)
    filepath = os.path.join(args.output_folder, filename)
    property.write_property_vector(filepath)
    """
    Generate all output files in case of a single sample size as input
    """
    if not isinstance(args.size, list):
        log.info(
            "{:15s} Sample percentage size: {:4.5s}".format("INPUT", str(args.size))
        )
        write_output_files(
            output_folder=args.output_folder,
            file_prefix=args.file_prefix,
            p_prop=property,
            s_prop=sampled_property,
            p_data=p_data,
            p=args.size,
            machine_learning=args.machine_learning,
        )
        plot_property(
            output_folder=args.output_folder,
            file_prefix=args.file_prefix,
            p_prop=property,
            s_prop=sampled_property,
            p=args.size,
        )
        log.info(
            "{:15s} Output files for selected sample size were generated successfully".format(
                "OUTPUT"
            )
        )
    """
    Create file with data report that includes important statistics about trajectory
    """
    filename = "{}_{}.json".format(args.file_prefix, "stats_report")
    filepath = os.path.join(args.output_folder, filename)
    p_data.property_data_report(filepath)
    log.info("{:15s} Statistics report was generated successfully".format("OUTPUT"))


def main(arg_list):
    sampling_workflow(arg_list)


def run():
    arg_list = sys.argv[1:]
    main(arg_list)


if __name__ == "__main__":
    run()
