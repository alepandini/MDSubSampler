from numpy import require
import mdss_geometrical_property
import mdss_pca_property
import mdss_sampler
import mdss_dissimilarity
import argparse
import sys
import operator as op

PROPERTY_CLASSES = [
    mdss_geometrical_property.RMSDProperty,
    mdss_geometrical_property.DistanceBetweenAtoms,
    mdss_geometrical_property.RadiusOfGyrationProperty,
    mdss_pca_property.TrjPCA,
    mdss_geometrical_property.DihedralAngles,
    mdss_geometrical_property.Angles,
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

DISSIMILARITY_CLASSES = [
    mdss_dissimilarity.Dissimilarity,
    mdss_dissimilarity.BhattaCoefficient,
    mdss_dissimilarity.KLDivergence,
    mdss_dissimilarity.PearsonCoefficient,
]
DISSIMILARITY_CLASS_MAPPING = dict(
    (dissimilarity.__name__, dissimilarity) for dissimilarity in DISSIMILARITY_CLASSES
)

# Function that creates the parser
def parse_args(arg_list):
    parser = argparse.ArgumentParser(description="Subsampler tool")
    parser.add_argument(
        "--traj",
        dest="trajectory_file",
        required=False,
        help="the path to the trajectory file",
    )
    parser.add_argument(
        "--top",
        dest="topology_file",
        required=False,
        help="the path to the topology file",
    )
    parser.add_argument(
        "--xvg",
        dest="xvg_file",
        required=False,
        help="the path to the xvg file containing the calculated property",
    )
    parser.add_argument(
        "--prefix",
        dest="file_prefix",
        required=True,
        help="the prefix for output files",
    )
    parser.add_argument(
        "--output-folder",
        dest="output_folder",
        required=True,
        help="the path to the output folder",
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
        "--strata-vector",
        dest="strata_vector",
        type=float,
        help=(
            "2D list of multiple layers and each layer is a set of labels "
            "for the frames according to the strata. This list should be in "
            "format of: [range(1, 20), range(20, 45),...] according to the "
            "frames in the protein trajectory"
        ),
    )
    parser.add_argument(
        "--strata-number",
        dest="strata_number",
        type=int,
        help=("The desired number of strata "),
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

    parser.add_argument(
        "--dissimilarity",
        dest="dissimilarity",
        choices=DISSIMILARITY_CLASS_MAPPING.keys(),
        help="Dissimilarity between distributions",
    )

    args = parser.parse_args(arg_list)

    if args.property in [
        mdss_geometrical_property.DistanceBetweenAtoms.__name__,
        mdss_geometrical_property.DihedralAngles.__name__,
        mdss_geometrical_property.Angles.__name__,
    ]:
        args.atom_selection = args.atom_selection.split(",")

    if args.xvg_file is not None and (
        args.trajectory_file is not None or args.topology_file is not None
    ):
        parser.print_help()
        print()
        print("Either provide an xvg file, or both trajectory and topology files")
        sys.exit(1)
    elif args.xvg_file is None:
        if args.trajectory_file is None or args.topology_file is None:
            parser.print_help()
            print()
            print("Both trajectory and topology files are required")
            sys.exit(1)

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
    require_sampler_argument(mdss_sampler.StratifiedSampler, "strata_vector")
    require_sampler_argument(mdss_sampler.UniformSampler, "strata_number")
    require_sampler_argument(mdss_sampler.BootstrappingSampler, "number_of_iterations")

    require_property_argument(
        mdss_geometrical_property.DistanceBetweenAtoms, "atom_selection"
    )
    require_property_argument(
        mdss_geometrical_property.DihedralAngles, "atom_selection"
    )
    require_property_argument(mdss_geometrical_property.Angles, "atom_selection")

    return args
