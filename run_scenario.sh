#!/usr/bin/env bash

python mdss_scenario.py --traj data/MD01_1lym_example_fit_short.xtc --top data/MD01_1lym_example.gro --prefix "$1"
