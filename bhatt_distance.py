

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
import seaborn as sns

# Code to read .XTC (trajectory) and .PDB(topology) file 
config = configparser.ConfigParser()
config.read("SAMPLE1.INI",)
path = config['PROTEINFILE']
trajectory_filename = path["BasePath"]+path["trajectory"]
topology_filename = path["BasePath"]+path["topology"]

u = mda.Universe(topology_filename,trajectory_filename, topology_format='PDB')

nterm = u.select_atoms('name N')[0]
cterm = u.select_atoms('name OC2')[-1]

bb = u.select_atoms('name CA')  # a selection (AtomGroup)

rgyr1 = []
time1 = []
d1 = []
r1 = []
for ts in u.trajectory[:5]:     # iterate through all frames
    r = cterm.position - nterm.position # end-to-end vector from atom positions
    d = numpy.linalg.norm(r)  # end-to-end distance
    time = u.trajectory.time
    rgyr = bb.radius_of_gyration()  # method of AtomGroup
    
    #print("frame = {0}: d = {1} A,r = {2}, Rgyr = {3} A, time = {4}".format(ts.frame, d, r, rgyr, time))
    time1.append(u.trajectory.time)
    rgyr1.append(bb.radius_of_gyration()) 
    d1.append(numpy.linalg.norm(r))
    r1.append(cterm.position - nterm.position)
    
rgyr_df = pd.DataFrame(rgyr1, columns=['Radius of gyration (A)'], index=time1)
d_df = pd.DataFrame(d1, columns=['Distance'], index=time1)
r_df = pd.DataFrame(r1, columns=['Vector1','Vector2','Vector3'], index=time1)
rgyr_df.index.name = 'Time (ps)'

frames = [d_df,rgyr_df,r_df]

result = pd.concat(frames, join='inner', axis=1)


h = result.to_numpy();

def mean( hist ):
    mean = 0.0;
    for i in hist:
        mean += i;
    mean/= len(hist);
    return mean;

def bhatta ( hist1,  hist2):
    # calculate mean of hist1
    h1_ = mean(hist1);

    # calculate mean of hist2
    h2_ = mean(hist2);

    # calculate score 
    score = 0;
    for i in range(len(h)):
        score += math.sqrt( hist1[i] * hist2[i] );
    # print h1_,h2_,score;
    score = math.sqrt( 1 - ( 1 / math.sqrt(h1_*h2_*8*8) ) * score );
    return score;

# generate and output scores
print("----------------------")
print("Bhattacharyya distance")
print("----------------------")
scores = [];
for i in range(len(h)):
    score = [];
    for j in range(len(h)):
        score.append( bhatta(h[i],h[j]) );
    scores.append(score);

for i in scores:
    print (i)
