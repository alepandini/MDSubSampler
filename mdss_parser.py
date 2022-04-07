from numpy import require
import mdss_protein_data
import mdss_property
import mdss_sampler
import argparse
import sys
import operator as op
import os

PROPERTY_CLASSES = [
    mdss_property.DistanceProperty,
    mdss_property.RadiusOfGyrationProperty,
    mdss_property.PCA,
    mdss_property.DihedralAngles,
    mdss_property.Angles,
]
PROPERTY_CLASS_MAPPING = dict(
    (property.__name__, property) for property in PROPERTY_CLASSES
)

SAMPLER_CLASSES = [
    mdss_sampler.RandomSampler,
    mdss_sampler.UniformSampler,
    mdss_sampler.StratifiedSampler,
    mdss_sampler.BootstrappingSampler,
]
SAMPLER_CLASS_MAPPING = dict((sampler.__name__, sampler) for sampler in SAMPLER_CLASSES)


# Function that creates the parser
def parse_args():
    parser = argparse.ArgumentParser(description="Subsampler tool")
    parser.add_argument(
        "--traj",
        dest="trajectory_file",
        required=True,
        help="the path to the trajectory file",
    )
    parser.add_argument(
        "--top",
        dest="topology_file",
        required=True,
        help="the path to the topology file",
    )
    parser.add_argument(
        "--prefix",
        dest="file_prefix",
        required=True,
        help="the prefix for output files",
    )
    parser.add_argument(
        "--output",
        dest="output_file",
        required=True,
        help="the path to the output file",
    )

    parser.add_argument(
        "--property",
        dest="property",
        choices=PROPERTY_CLASS_MAPPING.keys(),
        help="Property",
    )
    parser.add_argument(
        "--atom-selection",
        dest="atom_selection",
        help=(
            "Atom selection for calculation of any geometric property, "
            "If a selection is not specified then the default CA is used"
        ),
    )
    parser.add_argument(
        "--sampler",
        dest="sampler",
        choices=SAMPLER_CLASS_MAPPING.keys(),
        help="Sampler",
    )
    parser.add_argument(
        "--seed-number", dest="seed_number", type=int, help="Seed number"
    )
    parser.add_argument(
        "--low", dest="low", type=float, help="Lower boundary of the output interval"
    )
    parser.add_argument(
        "--high", dest="high", type=float, help="Higher boundary of the output interval"
    )
    parser.add_argument("--dtype", dest="dtype", type=int, help="dtype")
    parser.add_argument(
        "--layers",
        dest="layers",
        type=float,
        help=(
            "2D list of multiple layers and each layer is a set of labels "
            "for the frames according to the strata. This list should be in "
            "format of: [range(1, 20), range(20, 45),...] according to the "
            "frames in the protein trajectory"
        ),
    )

    parser.add_argument(
        "--n-iterations",
        dest="number_of_iterations",
        type=int,
        help="Number of iterations",
    )
    parser.add_argument(
        "--size",
        dest="size",
        type=int,
        help="Sample size",
    )

    args = parser.parse_args()

    if args.property in [
        mdss_property.DistanceProperty.__name__,
        mdss_property.DihedralAngles.__name__,
        mdss_property.Angles.__name__,
    ]:
        args.atom_selection = args.atom_selection.split(",")

    def require_sampler_argument(sampler, argument):
        value = op.attrgetter(argument)(args)
        if args.sampler == sampler.__name__ and value is None:
            parser.print_help()
            print()
            print("{} is required for {}".format(argument, sampler.__name__))
            sys.exit(1)

    def require_property_argument(property, argument):
        value = op.attrgetter(argument)(args)
        if args.property == property.__name__ and value is None:
            parser.print_help()
            print()
            print("{} is required for {}".format(argument, property.__name__))
            sys.exit(1)

    require_sampler_argument(mdss_sampler.RandomSampler, "seed_number")
    require_sampler_argument(mdss_sampler.UniformSampler, "low")
    require_sampler_argument(mdss_sampler.UniformSampler, "high")
    require_sampler_argument(mdss_sampler.StratifiedSampler, "layers")
    require_sampler_argument(mdss_sampler.BootstrappingSampler, "number_of_iterations")

    require_property_argument(mdss_property.DistanceProperty, "atom_selection")
    require_property_argument(mdss_property.DihedralAngles, "atom_selection")
    require_property_argument(mdss_property.Angles, "atom_selection")

    return args


if __name__ == "__main__":
    args = parse_args()
    print(args)

    file_prefix = args.file_prefix
    output = args.output
    if not os.path.exists(args.output):
        os.mkdirs(args.output)

    p_data = mdss_protein_data.ProteinData(
        args.trajectory_file, args.topology_file, config_parameters=None
    )

    ##################
    property_class = PROPERTY_CLASS_MAPPING[args.property]
    property = property_class(p_data, args.atom_selection)
    sampler_class = SAMPLER_CLASS_MAPPING[args.sampler]
    if args.sampler == "UniformSampler":
        sampler = sampler_class(p_data, args.low, args.high, args.dtype)
    elif args.sampler == "StratifiedSampler":
        sampler = sampler_class(p_data, args.layers)
    elif args.sampler == "BootstrappingSampler":
        sampler = sampler_class(p_data, args.number_of_iterations)
    ##################

    # Method that uses the user input for property calculation, sampling method
    # and size of sample and returns a subsample trajectory along with a log
    # file with diagnostics for the particular property.
    def run_subsampler(protein_data, property_class, sampler_class):
        print(
            f"Calculating {property_class.display_name}"
        )  # Will be replaced with log file
        print(f"Applying {property_class.display_name}")
        property.calculate_property()  # add this to the distribution class and just call the class here
