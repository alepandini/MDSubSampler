#!/usr/bin/env python
from mdss import main
import sys

if __name__ == "__main__":
    args = sys.argv[1:]
    trj_filename = args[0]
    top_filename = args[1]
    out_prefix = args[2]

    out_dir = "testing"
    property = "RMSDProperty"
    selection = "name CA"
    sampler = "RandomSampler"
    size = "100"
    dissimilarity = "Dissimilarity"

    size_list = [100, 200, 300, 400, 500, 600]
    for size in size_list:

        dissimilarity_score = main(
            [
                "--traj",
                trj_filename,
                "--top",
                top_filename,
                "--prefix",
                out_prefix,
                "--output-folder",
                out_dir,
                "--property",
                property,
                "--atom-selection",
                selection,
                "--sampler",
                sampler,
                "--seed-number",
                "1999",
                "--size",
                str(size),
                "--dissimilarity",
                dissimilarity,
            ]
        )
        print(dissimilarity_score)
        if dissimilarity_score < 0.001:
            break


# python scenario.py data/user.xtc data/user.gro testing
