#!/usr/bin/env python
from ctypes.wintypes import SIZE
from mdss import sampling_workflow
import src.mdss.utilities as u
import sys

OUT_DIR = "testing"
PROPERTY = "RMSD"
SELECTION = "name CA"
SAMPLER = "RandomSampler"
SIZE = "100"
DISSIMILARITY = "Dissimilarity"


def percentage(part, whole):
    percentage = 100 * float(part) / float(whole)
    return percentage


def main(trj_filename, top_filename, out_prefix):

    traj_size = u.check_trajectory_size(trj_filename, top_filename)
    perc_90 = percentage(90, traj_size)
    perc_80 = percentage(80, traj_size)
    perc_70 = percentage(70, traj_size)
    perc_60 = percentage(60, traj_size)
    perc_50 = percentage(50, traj_size)
    perc_40 = percentage(40, traj_size)

    size_list = [perc_90, perc_80, perc_70, perc_60, perc_50, perc_40]
    for size in size_list:

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
                "--seed-number",
                "1999",
                "--size",
                str(SIZE),
                "--dissimilarity",
                DISSIMILARITY,
            ]
        )
        print("for {} % sample of full trajectory:".format(size))
        print(dissimilarity_score)
        if dissimilarity_score < 0.5:
            break


if __name__ == "__main__":
    args = sys.argv[1:]
    trj_filename = args[0]
    top_filename = args[1]
    out_prefix = args[2]
    main(trj_filename, top_filename, out_prefix)

# python scenario.py data/input.xtc data/input.gro testing
