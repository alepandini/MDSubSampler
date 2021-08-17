import MDAnalysis as mda
import MDAnalysis
from MDAnalysis.analysis import rms
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import distances
from MDAnalysis.tests.datafiles import PSF, DCD
import numpy.linalg
import numpy as np
import numpy as lyzo
import pandas as pd
import re
import csv
import os
import pdb
import math
from Bio import SeqIO
#from pudb import set_trace
from Bio.PDB import *
from glob import glob
import configparser
from sklearn import preprocessing
import matplotlib.pyplot as plt
import pdb;
import random
from pandas.plotting import scatter_matrix
# from bhatta_dist import bhatta_dist

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
        #print(cAlphaSample) 

    def _distance_measure(self,u):
        # can access via atom name
        # we take the first atom named N and the last atom named C
        nterm = u.select_atoms('name N')[0]
        cterm = u.select_atoms('name OC2')[-1]

        bb = u.select_atoms('name CA')  # a selection (AtomGroup)
        # printing coordinates , based on that information how to pick the frame ?
        print("----------------------")
        print("COORDINATES INFORMATION")
        print("----------------------")
        print(bb.positions)
        rgyr1 = []
        time1 = []
        d1 = []
        r1 = []
        print("----------------------")
        print("TRAJECTORY INFORMATION")
        print("----------------------")
        for ts in u.trajectory[:6]:     # iterate through all frames
            r = cterm.position - nterm.position # end-to-end vector from atom positions
            d = numpy.linalg.norm(r)  # end-to-end distance
            time = u.trajectory.time
            rgyr = bb.radius_of_gyration()  # method of AtomGroup
            coordinates = np.random.rand(len(u.atoms), 3)
            print("frame = {0}: d = {1} A,r = {2}, Rgyr = {3} A, time = {4}".format(ts.frame, d, r, rgyr, time))
            time1.append(u.trajectory.time)
            rgyr1.append(bb.radius_of_gyration()) 
            d1.append(numpy.linalg.norm(r))
            r1.append(cterm.position - nterm.position) 
        # saving each attibutes in dataframe  
        rgyr_df = pd.DataFrame(rgyr1, columns=['Radius of gyration (A)'], index=time1)
        d_df = pd.DataFrame(d1, columns=['Distance'], index=time1)
        r_df = pd.DataFrame(r1, columns=['Vector1','Vector2','Vector3'], index=time1)
        rgyr_df.index.name = 'Time (ps)'
        #  concat all dataframe into one to create single input file        
        frames = [d_df,rgyr_df,r_df]
        result = pd.concat(frames, join='inner', axis=1)
        h = result.to_numpy();
        print("----------------------")
        print("DISTANCE MEASURE")
        print("----------------------")
        print(h)
        # save dataframe as CSV file
        #result.to_csv('file1.csv')
        self._statistical_analysis(result)
       
        # calling method to perfrom sampling 
        self._trajectroy_sampling(result)
        # self._bhatt_distance(h)
        # Ploting graph for radius of gyration against time 
        # rgyr_df.plot(title='Radius of gyration')
        # plt.savefig("rgyrgraph.jpg")
            
    def _statistical_analysis(self,result):
        print("----------------------")
        print("STATISTICAL ANALYSIS")
        print("----------------------")
        print(result.describe())
        # histogram for full trajectory
        result.hist()
        plt.savefig("FullTraj_Hist.jpg")
        # commented due to dynamic image path giving error 
        #plt.savefig(path["ImagePath"]+"FullTraj_Hist.jpg")
        # scatter matrix for full trajectory
        scatter_matrix(result, alpha=0.2, figsize=(6, 6), diagonal='kde')
        plt.savefig("FullTraj_Scatter.jpg")
        # commented due to dynamic image path giving error 
        # plt.savefig(path["ImagePath"]+"FullTraj_Scatter.jpg")
    
    def _trajectroy_sampling(self,result):
        print("----------------------")
        print("SAMPLING")
        print("----------------------")
        # random sampling 
        # axis=1  randomly sample columns and axis=0 for rows
        # seed for the random number generator can be specified using randon_state , same rows & columns retuns each time
        subset = result.sample(n=5,random_state=0)
        print("random sampling",subset)
        # histogram for subset 
        subset.hist()
        # commented due to dynamic image path giving error 
        # plt.savefig(path["ImagePath"]+"Sample_Hist.jpg")
        plt.savefig("Sample_Hist.jpg")
        # scatter matrix for subset
        scatter_matrix(subset, alpha=0.2, figsize=(6, 6), diagonal='kde')
        # commented due to dynamic image path giving error 
        plt.savefig("Sample_Scatter.jpg")
        #plt.savefig(path["ImagePath"]+"Sample_Scatter.jpg")
        
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
        # convert dataframe to numpy array
        
        #print(coordinates)
       
        print(os.getcwd())
        
if __name__=='__main__':
    main()
