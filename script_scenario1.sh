#!/usr/bin/env bash
"""
Drug design scenario: 

Problem: Select protein structures from trajectory for drug design. 
Goal: Extract subset of frames with similar distribution of pocket size.

Sampling method: Random Sampling
Property Calculation: RMSD
"""

# python mdss_parser.py \
#   --traj "data/user.xtc" \
#   --top "data/user.gro" \
#   --prefix "user" \
#   --output-folder "data/results" \
#   --property='RMSDProperty' \
#   --atom-selection='name CA' \
#   --sampler='RandomSampler' \
#   --seed-number=1999 \
#   --size=100 \
#   --distance='Distance'


#Testing for git