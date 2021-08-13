#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from MDAnalysis.tests.datafiles import PSF, DCD   # test trajectory
import numpy.linalg

import numpy
import math
from sklearn import preprocessing
import numpy as np
import pandas as pd
import csv
import re
import configparser
import MDAnalysis as mda
import MDAnalysis
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt

# CODE TO USE LYZO .XTC AND .GRO FILE 
#config = configparser.ConfigParser()
#config.read("SAMPLE1.INI",)
#path = config['PROTEINFILE']
#trajectory_filename = path["BasePath"]+path["trajectory"]
#topology_filename = path["BasePath"]+path["topology"]
# u = mda.Universe(topology_filename,trajectory_filename, topology_format='GRO')
#CODE ENDS HERE

u = MDAnalysis.Universe(PSF,DCD)  # always start with a Universe

# can access via segid (4AKE) and atom name
# we take the first atom named N and the last atom named C

nterm = u.select_atoms('segid 4AKE and name N')[0]
cterm = u.select_atoms('segid 4AKE and name C')[-1]
#done some ammendement according to PDB file we have but if we use beloe code - gives below error 
#IndexError: index 0 is out of bounds for axis 0 with size 0
#nterm = u.select_atoms('segid 1MET and name N')[0]
#cterm = u.select_atoms('segid 164LEU and name C')[-1]

bb = u.select_atoms('name CA')  # a selection (AtomGroup)

for ts in u.trajectory:     # iterate through all frames
    r = cterm.position - nterm.position # end-to-end vector from atom positions
    d = numpy.linalg.norm(r)  # end-to-end distance
    rgyr = bb.radius_of_gyration()  # method of AtomGroup
    print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
          ts.frame, d, rgyr))

    print("end-to-end vector::",r)
    print("end-to-end distance::",d)