#!/usr/bin/env python
from mdss import main, sampling_workflow
from mdss_utilities import check_trajectory_size as tr
import sys

out_dir = "testing"
property = "RMSDProperty"
selection = "name CA"
sampler = "RandomSampler"
size = "100"
dissimilarity = "Dissimilarity"


def percentage(part, whole):
    percentage = 100 * float(part) / float(whole)
    return percentage


def main(trj_filename, top_filename, out_prefix):

    traj_size = tr.check_trajectory_size(trj_filename, top_filename)
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
        print("for {} percent sample of full trajectory:".format(size))
        print(dissimilarity_score)
        if dissimilarity_score < 0.001:
            break


if __name__ == "__main__":
    arg_list = sys.argv[1:]
    main(arg_list)
    trj_filename = args[0]
    top_filename = args[1]
    out_prefix = args[2]

# python scenario.py data/user.xtc data/user.gro testing
