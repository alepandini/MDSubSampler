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
python3.10 mdss_parser.py \
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

