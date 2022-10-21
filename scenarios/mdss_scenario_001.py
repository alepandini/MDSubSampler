#!/usr/bin/env python
from mdss import sampling_workflow
import sys

OUT_DIR = "testing"
PROPERTY = "RMSDProperty"
SELECTION = "name CA"
SAMPLER = "RandomSampler"
SIZE = "100"
DISSIMILARITY = "Dissimilarity"
FIT = "False"


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
            "--fit",
            FIT,
            "--atom-selection",
            SELECTION,
            "--sampler",
            SAMPLER,
            "--seed-number",
            "1999",
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
