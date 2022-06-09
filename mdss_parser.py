from numpy import require
import mdss_protein_data
import mdss_property
import mdss_sampler
import mdss_distance
import argparse
import sys
import operator as op
import os
import mdss_distribution

PROPERTY_CLASSES = [
    mdss_property.RMSDProperty,
    mdss_property.DistanceBetweenAtoms,
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

DISTANCE_CLASSES = [
    mdss_distance.Distance,
    mdss_distance.BhattaDistance,
    mdss_distance.KLDiverDistance,
    mdss_distance.PearsonDictDistance,
]
DISTANCE_CLASS_MAPPING = dict(
    (distance.__name__, distance) for distance in DISTANCE_CLASSES
)

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

    parser.add_argument(
        "--distance",
        dest="distance",
        choices=DISTANCE_CLASS_MAPPING.keys(),
        help="Distance between distributions",
    )

    args = parser.parse_args()

    if args.property in [
        mdss_property.DistanceBetweenAtoms.__name__,
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

    require_property_argument(mdss_property.DistanceBetweenAtoms, "atom_selection")
    require_property_argument(mdss_property.DihedralAngles, "atom_selection")
    require_property_argument(mdss_property.Angles, "atom_selection")

    return args


####################################

if __name__ == "__main__":
    args = parse_args()
    print(args)

    file_prefix = args.file_prefix
    output_folder = args.output_folder
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    p_data = mdss_protein_data.ProteinData(
        args.trajectory_file, args.topology_file, config_parameters=None
    )

    property_class = PROPERTY_CLASS_MAPPING[args.property]
    property = property_class(p_data, args.atom_selection)
    property.calculate_property()

    sampler_class = SAMPLER_CLASS_MAPPING[args.sampler]
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

    distance_class = DISTANCE_CLASS_MAPPING[args.distance]
    distance = distance_class(property, property_sample)
    print("Distance: {}".format(distance.calculate_distance()))
    distribution = mdss_distribution.DistributionDistanceSimple(
        property, property_sample, distance
    )
    distribution.simple_distance_between_distributions()

    distrib = mdss_distribution.DistributionDistance(
        property, property_sample, distance
    )
    filename = "{}_{}_{}.dat".format(
        file_prefix, property_class.display_name, distance_class.display_name
    )
    filepath = os.path.join(args.output_folder, filename)

    property.write_property_vector(filepath)

    filename_sample = "{}_{}_sample_{}.dat".format(
        file_prefix, property_class.display_name, distance_class.display_name
    )
    filepath_sample = os.path.join(args.output_folder, filename_sample)

    property_sample.write_property_vector(filepath_sample)

    filename_plot = "{}_{}_{}_plot.png".format(
        file_prefix, property_class.display_name, distance_class.display_name
    )
    filepath_plot = os.path.join(args.output_folder, filename_plot)
    property.plot_property(filepath, filepath_plot)

    filename_sample_plot = "{}_{}_sample_{}_plot.png".format(
        file_prefix, property_class.display_name, distance_class.display_name
    )
    filepath_sample_plot = os.path.join(args.output_folder, filename_sample_plot)
    property_sample.plot_property(filepath_sample, filepath_sample_plot)
    ####################################

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
        property.calculate_property()  # add this to the distribution class and just call the class here


# python mdss_parser.py --traj "data/MD01_1lym_example_fit_short.xtc" --top "data/MD01_1lym_example.gro" --prefix "001" --output-folder "data/results" --property='RMSDProperty' --atom-selection='name CA' --sampler='RandomSampler' --seed-number=1999 --size=100 --distance='Distance'
