import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy.linalg
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
# This class designed to read trajectory file , topology file and calculate distace measure for particular ATOM (CA).Code is aim # to adapt user atom selection
class Protein_Data:
    
    def __init__(self,trajectory_filename,topology_filename):
       # pdb.set_trace()
        self.trajectory_data = self._read_trajectory(trajectory_filename,topology_filename)
        self.transform_data = self._transform_trajectory()
        self.select_atom = self._select_atoms()
        self.normalise_traj = self._normalise_trajectory(topology_filename)
        
        print("4 constructor exit")
        # This method take trajectory and topology file as input and return Universe object which can 
        # be used for further calculation 
    def _read_trajectory(self,trajectory_filename,topology_filename):
        print("2 read")
        u = mda.Universe(topology_filename,trajectory_filename,permissive=False, topology_format='PDB')
        # Call distance measure function and taking universe object as parameter to provide end-to-end vector 
        # and end-to-end vector distance
        self._distance_measure(u)
        return u
    def _transform_trajectory(self):
        print("3 transform")
         # This function currently normalise topology file using static file path. 
    def _normalise_trajectory(self,topology_filename):
        
        traj_filename = "/Users/riktapatel/Documents/Disertation/Data/MD01_1lym_example.gro"
         
         #  skiprows and skipfooter parameter used to remove top and bottom rows from topology file 
         #  delim_whitespace specifies whether or not whitespace ,no delimiter parameter if option is true.
         #  encoding to use for UTF when reading/writing files
         #  List of column names to use. if file contains a header row, then header=0 to override column name      
       
        Cov = pd.read_table(traj_filename, 
                    header=0, 
                    skiprows=2,
                    skipfooter=1,
                    encoding='utf-8',
                    delim_whitespace=True,
                    names=["A", "B", "C", "D", "E", "F", "G","H","I"]
                    )
        # axis=1 when drops column 
        Dropresult = Cov.drop(['A', 'B'], axis=1)
        # below lines of code used to normalised topology data 
        d = preprocessing.normalize(Dropresult, axis=0)
        scaled_df = pd.DataFrame(d, columns=["C", "D", "E", "F", "G","H","I"])
        scaled_df.head()
        cAlphaSeq = np.array(scaled_df)
        # below lines of code used to sample topology data
        cAlphaSample = scaled_df.sample(frac=0.5, replace=True, random_state=1) 
        print(cAlphaSample) 

    def _distance_measure(self,u):
        # can access via atom name
        # we take the first atom named N and the last atom named C
        nterm = u.select_atoms('name N')[0]
        cterm = u.select_atoms('name OC2')[-1]

        bb = u.select_atoms('name CA')  # a selection (AtomGroup)

        for ts in u.trajectory:     # iterate through all frames
            r = cterm.position - nterm.position # end-to-end vector from atom positions
            d = numpy.linalg.norm(r)  # end-to-end distance
            rgyr = bb.radius_of_gyration()  # method of AtomGroup
            print("frame = {0}: d = {1} A, Rgyr = {2} A".format(ts.frame, d, rgyr))

            print("end-to-end vector::",r)
            print("end-to-end distance::",d)
    def _select_atoms(self):
        print("4 select atoms")
        
def main():
        #initialize the parser 
        config = configparser.ConfigParser()
        # load the configuration file
        config.read("SAMPLE1.INI",)
        # read values from a relevent section header
        path = config['PROTEINFILE']
        trajectory_filename = path["BasePath"]+path["trajectory"]
        topology_filename = path["BasePath"]+path["topology"]
        prodata = Protein_Data(trajectory_filename,topology_filename)
        traj = prodata.trajectory_data
        # selecting atom called "C Aplha" from trajectory file
        Calpha = traj.select_atoms("name CA")
        len(traj.atoms)
        # printing coordinated for C Alpha
        coordinates = np.random.rand(len(traj.atoms), 3)
        print(coordinates)
        print(os.getcwd())
if __name__=='__main__':
    main()
