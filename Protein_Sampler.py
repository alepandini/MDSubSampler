import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
import numpy as lyzo
import pandas as pd
import re
import csv
import os
import pdb
from Bio import SeqIO
#from pudb import set_trace
from Bio.PDB import *
from glob import glob
import configparser
from sklearn import preprocessing

class Protein_Data:
    def __init__(self,trajectory_filename,topology_filename):
       # pdb.set_trace()
        self.trajectory_data = self._read_trajectory(trajectory_filename,topology_filename)
        self.transform_data = self._transform_trajectory()
        self.select_atom = self._select_atoms()
        self.normalise_traj = self._normalise_trajectory(trajectory_filename)
        print("4 constructor exit")
    def _read_trajectory(self,trajectory_filename,topology_filename):
        print("2 read")
        u = mda.Universe(topology_filename,trajectory_filename,permissive=False, topology_format='GRO')
        return u
    def _transform_trajectory(self):
        print("3 transform")
    def _normalise_trajectory(self,trajectory_filename):
        traj_filename = "/Users/riktapatel/Documents/Disertation/Data/MD01_1lym_example_norm.gro"
        Cov = pd.read_table( traj_filename, 
                    header=0, 
                    skiprows=2,
                    skipfooter=1,
                    encoding='utf-8',
                    delim_whitespace=True,
                    names=["A", "B", "C", "D", "E", "F", "G","H","I"]
                    )
        print(Cov)
        Dropresult = Cov.drop(['A', 'B'], axis=1)
        print("testing",Dropresult)
        d = preprocessing.normalize(Dropresult, axis=0)
        scaled_df = pd.DataFrame(d, columns=["C", "D", "E", "F", "G","H","I"])
        scaled_df.head()
        cAlphaSeq = np.array(scaled_df)
        print("Printing sample size here")
        cAlphaSample = scaled_df.sample(frac=0.5, replace=True, random_state=1) 
        print(cAlphaSample) 
   # def print_coordinates(neighborList):
        #for y in neighborList:
        #print("6 th " , y.get_coord())
    def _select_atoms(self):
        print("4 select atoms")
        
def main():
    config = configparser.ConfigParser()
    config.read("SAMPLE1.INI",)
    path = config['PROTEINFILE']
    trajectory_filename = path["BasePath"]+path["trajectory"]
    topology_filename = topology_filename = path["BasePath"]+path["topology"]
    #print("1 main")
    prodata = Protein_Data(trajectory_filename,topology_filename)
    #print(pdata1)
    #print("object",pdata1)
    traj = prodata.trajectory_data
     #normdata = traj.normalise_trajectory
    #print("testing",traj)
    Calpha = traj.select_atoms("name CA")
    print(Calpha)
    len(traj.atoms)
    coordinates = np.random.rand(len(traj.atoms), 3)
    print(coordinates)
    #lyzo_atoms = lyzo.atoms.positions
    #print("lyzo atoms",lyzo_atoms)
   # calpha = lyzo.select_atoms("name CA")
     # lyzo_Cord = np.random.rand(len(lyzo.atoms),3)
     # print(lyzo_Cord)
     # lyzo.load_new(lyzo_Cord,format=MemoryReader)
     # print("atom positions",lyzo.atoms.positions)
    
    print(os.getcwd())
    #glob.glob('subdir/*.pdb')

if __name__=='__main__':
    main()