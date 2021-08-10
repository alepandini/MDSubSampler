#!/usr/bin/env python
# coding: utf-8

# In[8]:


import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
import numpy as lyzo
import os
import pdb
#from Bio import SeqIO
import configparser
#from pudb import set_trace;
from Bio.PDB import *
from glob import glob

class Fit_Trajectory:
    
    def __init__(self,trajectory_filename,topology_filename):
       
        self.fit_traj = self._fit_traj()
        print("4 constructor exit")
        
    def fit_traj():
        config = configparser.ConfigParser()
        config.read("SAMPLE1.INI",)
        path = config['PROTEINFILE']
        print("trajectoryfile :",path["BasePath"]+path["trajectory"])
        trajectory_filename = path["BasePath"]+path["trajectory"]
        topology_filename = path["BasePath"]+path["topology"]
        traj_fit_filename = path["BasePath"]+path["rmsfitfile"]
   
        print(trajectory_filename)
        print(topology_filename)
        LyzoRef = mda.Universe(topology_filename)
        LyzoData = mda.Universe(topology_filename, trajectory_filename,permissive=False)
        trajAlign = align.AlignTraj(LyzoData, LyzoRef, filename = traj_fit_filename)
        alignment.run()
    
    def main():
        
         print("1 main")
            
    if __name__ == '__main__':
        main()
        fit_traj()


# In[ ]:




