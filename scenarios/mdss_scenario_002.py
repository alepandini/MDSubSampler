#!/usr/bin/env python
from mdss import sampling_workflow
import sys

OUT_DIR = "testing"
PROPERTY = "RMSDProperty"
SELECTION = "resid 55:127"
SAMPLER = "UniformSampler"
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


if __name__ == "__main__":
    args = sys.argv[1:]
    trj_filename = args[0]
    top_filename = args[1]
    out_prefix = args[2]
    main(trj_filename, top_filename, out_prefix)
