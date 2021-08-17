#!/usr/bin/env python
# coding: utf-8

# In[8]:


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
from pandas.plotting import scatter_matrix

# Code to read .XTC (trajectory) and .PDB(topology) file 
config = configparser.ConfigParser()
config.read("SAMPLE1.INI",)
path = config['PROTEINFILE']
trajectory_filename = path["BasePath"]+path["trajectory"]
topology_filename = path["BasePath"]+path["topology"]
rmsf_filename = path["BasePath"]+path["rmsffitfile"]
rmsd_filename = path["BasePath"]+path["rmsdfitfile"]


u = mda.Universe(topology_filename,trajectory_filename, topology_format='PDB')
# RMSF CODE 
resid, rmsf = numpy.loadtxt(rmsf_filename, unpack=True)

fig = plt.figure(figsize=(5,2.5))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)

ax.fill_between(resid, rmsf, color="red", linestyle="-", alpha=0.1)
ax.plot(resid, rmsf, color="red", linestyle="-")

# RMSD CODE 
t,rmsd = numpy.loadtxt(rmsd_filename, unpack=True)

fig = plt.figure(figsize=(5,2.5))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)

ax.fill_between(t,rmsd, color="blue", linestyle="-", alpha=0.1)
ax.plot(t,rmsd, color="blue", linestyle="-")


# In[ ]:





# In[ ]:




