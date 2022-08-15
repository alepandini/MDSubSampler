#!/bin/bash
if [ $# -ne 3 ] 
then
    echo "Incorrect number of arguments..."
    echo "Usage: mdss_scenario_001.sh <traj_filename> <top_filename> <output_prefix>"
    exit 1
fi

######################################################################
#   user defined parameters
TRJ_FILENAME=$1
TOP_FILENAME=$2
OUT_PREFIX=$3

######################################################################
#   preset parameters
OUT_DIR=$PWD
PROPERTY="RMSDProperty"
SELECTION="name CA"
SAMPLER="RandomSampler"
SIZE=100
DISTANCE="Distance"

######################################################################
#   wrapper to MDSS python call
python3.10 mdss.py \
  --traj $TRJ_FILENAME \
  --top $TOP_FILENAME \
  --prefix $OUT_PREFIX \
  --output-folder $OUT_DIR \
  --property=$PROPERTY \
  --atom-selection="$SELECTION" \
  --sampler=$SAMPLER\
  --seed-number=1999 \
  --size=$SIZE \
  --distance=$DISTANCE



# #!/usr/bin/env bash
# """
# Drug design scenario: 

# Problem: Select protein structures from trajectory for drug design. 
# Goal: Extract subset of frames with similar distribution of pocket size.

# Sampling method: Random Sampling
# Property Calculation: RMSD
# """

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
#   --distance='Dissimilarity'