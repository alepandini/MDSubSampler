"""
Drug design scenario: 

Problem: Select protein structures from trajectory for drug design. 
Goal: Extract subset of frames with similar distribution of pocket size.

Sampling method: Random Sampling
Property Calculation: RMSD
"""

import mdss_protein_data
import mdss_property
import mdss_sampler
import mdss_distance
import mdss_distribution
import argparse
import mdss_parser as p
import os

# Create output folder for results
args = p.parse_args()
file_prefix = args.file_prefix
output_folder = args.output_folder
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

# Import data - protein trajectory
p_data = mdss_protein_data.ProteinData(
    args.trajectory_file, args.topology_file, config_parameters=None
)

# Calculate the property - RMSD - of the trajectory frames
property_class = p.PROPERTY_CLASS_MAPPING[args.property]
property = property_class(p_data, args.atom_selection)
property.calculate_property()

# Subsample the protein trajectory with the desired sample size
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

# Calculate distance between full trajectory and sample trajectory
distance_class = p.DISTANCE_CLASS_MAPPING[args.distance]
distance = distance_class(property, property_sample)
print("Distance: {}".format(distance.calculate_distance()))
distribution = mdss_distribution.DistributionDistanceSimple(
    property, property_sample, distance
)
distribution.simple_distance_between_distributions()
distrib = mdss_distribution.DistributionDistance(
    property, property_sample, distance
)

# Output files and plots of full and sample trajectory
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
