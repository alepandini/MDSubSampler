"""

Scenario 003

Purpose:             Given a single MD trajectory sampling different conformational 
                     motions of the protein, select a subset of frames proportionally 
                     to the most frequent values of RMSD  

User:                Molecular dynamics user with basic understanding of coding 
    
Input:               Molecular dynamics trajectory 
                     Subsample size 

Sampling strategy:   Weighted random sampling 

Type of property:    Numerical continuous 

Property dependency: Dependent on single frame 

Criterion:           Dissimilarity between distributions of values of original 
                     and subsampled property 

Scenario:	
                    1. Read input trajectory and topology files 
                    2. Read subsample size 
                    3. Calculate RMSD distribution for input trajectory 
                    4. Weighted random sampling 
                    5. Calculate RMSD distribution for subsampled trajectory 
                    6. Calculate dissimilarity measure between distributions 

"""

#!/usr/bin/env python
from mdss.run import sampling_workflow
import sys

OUT_DIR = "testing"
PROPERTY = "RMSDProperty"
SELECTION = "name CA"
SAMPLER = "WeightedSampler"
STRATA_NUMBER = "200"
SIZE = "100"
DISSIMILARITY = "Bhattacharya"


def main(trj_filename, top_filename, out_prefix):

    dissimilarity_score = sampling_workflow(
        [
            "--traj",
            trj_filename,
            "--top",
            top_filename,
            "--prefix",
            out_prefix,
            "--output-folder",
            OUT_DIR,
            "--property",
            PROPERTY,
            "--atom-selection",
            SELECTION,
            "--sampler",
            SAMPLER,
            "--strata-number",
            STRATA_NUMBER,
            "--size",
            str(SIZE),
            "--dissimilarity",
            DISSIMILARITY,
        ]
    )
    print(dissimilarity_score)


def run():
    args = sys.argv[1:]
    trj_filename = args[0]
    top_filename = args[1]
    out_prefix = args[2]
    main(trj_filename, top_filename, out_prefix)


if __name__ == "__main__":
    run()
