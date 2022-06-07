"""
Drug design scenario: 

Select structures for drug design. 
Goal: Extract subset of frames with similar distribution of pocket size.

Sampling method: Random Sampling
Property Calculation: RMSD
"""

import mdss_protein_data
import mdss_property
import mdss_sampler
import mdss_distance
import mdss_distribution
import os
import mdss_parser as p
from mdss_parser import argparse

# Creating an output folder for the results
# and loading the data
args = p.parse_args()
print(args)

file_prefix = args.file_prefix
output_folder = args.output_folder
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

p_data = mdss_protein_data.ProteinData(
    args.trajectory_file, args.topology_file, config_parameters=None
)

# Subsample the protein trajectory
property_class = p.PROPERTY_CLASS_MAPPING[args.property]
property = property_class(p_data, args.atom_selection)
property_sample = property_class(p_data, args.atom_selection)
sampler_class = p.SAMPLER_CLASS_MAPPING[args.sampler]
