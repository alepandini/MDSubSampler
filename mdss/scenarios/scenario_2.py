"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Namir Oues
    
    This file is part of the MDSubSampler software 
    (https://github.com/alepandini/MDSubSampler).
    Copyright (c) 2023 Namir Oues and Alessandro Pandini.

    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>.


Scenario 002

Purpose:             Given a single MD trajectory sampling different conformations 
                     of a binding pocket, select a subset of frames with a range of 
                     geometries for pocket opening.

User:                Molecular dynamics user with basic understanding of coding 
    
Input:               Molecular dynamics trajectory 
                     List of atoms for RMSD calculation 
                     Reference structure of open (or close) state 
                     List of atoms for RMSD fit [optional] 
                     Subsample size and number of intervals along RMSD values 

Sampling strategy:   Uniform random sampling 

Type of property:    Numerical continuous 

Property dependency: Dependent on single frame 

Criterion:           Dissimilarity between distributions of values of original 
                     and subsampled property 

Scenario:	
                    1. Read input trajectory and topology files 
                    2. Read atoms lists, reference structure, subsample size and number
                       of intervals 
                    3. Calculate RMSD distribution for input trajectory 
                    4. Uniform random subsample 
                    5. Calculate RMSD distribution for subsampled trajectory 
                    6. Calculate dissimilarity measure between distributions 

"""

#!/usr/bin/env python
from mdss.run import sampling_workflow
import sys

OUT_DIR = "results"
PROPERTY = "RMSD"
SELECTION = "resid 120:160"
SAMPLER = "UniformSampler"
STRATA_NUMBER = "200"
SIZE = "10"
DISSIMILARITY = "Bhattacharyya"


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
            "--fit",
            # "machine-learning",
        ]
    )


def run():
    args = sys.argv[1:]
    trj_filename = args[0]
    top_filename = args[1]
    out_prefix = args[2]
    main(trj_filename, top_filename, out_prefix)


if __name__ == "__main__":
    run()
