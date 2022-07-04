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

    main(
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
            size,
            "--dissimilarity",
            dissimilarity,
        ]
    )
