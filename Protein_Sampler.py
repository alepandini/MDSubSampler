import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)
import numpy.linalg
import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
from sklearn import preprocessing
import os
import configparser
import random
import dictances 
from dictances import bhattacharyya, bhattacharyya_coefficient
import itertools 
from math import log2
from math import sqrt
from typing import Dict
from distances_utils import sort



class Protein_Data:
    
    def __init__(self,trajectory_filename,topology_filename, config_parameters):
        self.config_par = config_parameters
        self.trajectory_filename = trajectory_filename
        self.topology_filename = topology_filename
        self.trajectory_data = self._read_trajectory(self.trajectory_filename,self.topology_filename)
        self.n_frames = self.trajectory_data.trajectory.n_frames
        self.ca_atom_group = self._select_CA_atoms()
        self.masses = self.ca_atom_group.masses
        self.total_mass = np.sum(self.masses)
        self.property_dict = {}

    def _read_trajectory(self,trajectory_filename,topology_filename):
        trajectory_data = mda.Universe(topology_filename,trajectory_filename,permissive=True, topology_format='PDB')
        return trajectory_data

    def _select_CA_atoms(self):
        ca_atom_group = self.trajectory_data.select_atoms("name CA")
        return ca_atom_group

    def output_trj_summary(self):
        print("----------------------")
        print("TRAJECTORY INFORMATION")
        print("----------------------")
        print("n frames = {0}\nn atoms = {1}\nn CA atoms = {2}".format(self.n_frames, self.trajectory_data.trajectory.n_atoms, self.ca_atom_group.n_atoms))

    def statistical_analysis(self):
        print("----------------------")
        print("STATISTICAL ANALYSIS")
        print("----------------------")
        print(result.describe())
    def add_property(self, property_vector, property_name, sample_label):
        if property_name in self.property_dict:
            self.property_dict[property_name][sample_label] = property_vector
        else:
            self.property_dict[property_name] = {}
            self.property_dict[property_name][sample_label] = property_vector
            

class Property_RMSD:
    def __init__(self, protein_data, frame_list, atom_selection = "name CA"):
        protein_data.trajectory_data.trajectory[0]
        ref_coordinates = protein_data.trajectory_data.select_atoms(atom_selection).positions.copy()
        
        self.rmsd = []
        for frame in frame_list:
            protein_data.trajectory_data.trajectory[frame]
            self.rmsd.append(
                rms.rmsd(
                    protein_data.trajectory_data.select_atoms(atom_selection).positions,
                    ref_coordinates
                ) 
            )
class Property_RadiusOfGyration:
       def __init__(self,protein_data, frame_list, atom_selection = "name CA", verbose=True):
       
        protein_data.trajectory_data.trajectory[0]
       
        self.Rgyr = [] 
        for frame in frame_list:
            protein_data.trajectory_data.trajectory[frame] 
            self.Rgyr.append((protein_data.trajectory_data.trajectory.time, 
                         protein_data.trajectory_data.select_atoms(atom_selection).radius_of_gyration()))
          
class Bhatta_Distance: 
        
        def __init__(self,protein_data,rmsd_vector,sampled_rmsd_vector, verbose=True):
            
        
            for size in [10,50,100]:
                lst_vector = listToDict(rmsd_vector)
             
                lst_sample_vector = listToDict(sampled_rmsd_vector)
                
                for i,size in enumerate(lst_sample_vector):
                    
                    bhatta_coeff = bhattacharyya_coefficient(lst_vector, lst_sample_vector)
                    
                    bhatta_dist = bhattacharyya(lst_vector, lst_sample_vector) 
                
            print(size+1,bhatta_dist)
            print(size+1,bhatta_coeff)
        
            
def listToDict(lst):
        op = { i : lst[i] for i in range(0, len(lst) ) }
        return op 


            
class Frame_Sampler:
        def __init__(self, frame_list, seed_number = 1999):
            random.seed(seed_number)
            self.frame_list = frame_list
            self.sampled_frame_list = None

        def sample(self, size):
            self.sampled_frame_list = random.sample(self.frame_list, size)

def get_config_parameters(config_filename):
        config = configparser.ConfigParser()
        config.read(config_filename,)
        config_par = config['PROTEINFILE']
        return config_par
def kbinDist(a):   
            data1=np.array(a)
            print(len(data1))
            data1 = data1.reshape((len(data1),1))
            print(data1)
            kbins = KBinsDiscretizer(n_bins=10, encode='ordinal', strategy='uniform')
            data_trans1 = kbins.fit_transform(data1)
            return data_trans1
def main():
        config_par = get_config_parameters("SAMPLE1.INI")
        
        trajectory_filename = config_par["BasePath"]+config_par["trajectory"]
        topology_filename = config_par["BasePath"]+config_par["topology"]
        #       creates a Protein_Data object
        pro_data = Protein_Data(trajectory_filename, topology_filename, config_par)

        pro_data.output_trj_summary()
        #       calculates an example property
        rmsd_vector = Property_RMSD(pro_data, range(pro_data.n_frames)).rmsd
        
        rgyr_vector = Property_RadiusOfGyration(pro_data, range(pro_data.n_frames))
        
        pro_data.add_property(rmsd_vector, "RMSD", "reference")
        pro_data.add_property(rgyr_vector, "Radius Of Gyration", "reference")
        #       creates a Frame_Sampler object
        frame_sampler = Frame_Sampler(range(pro_data.n_frames))
        #       for different values of sample size, the sampler randomly selects frames
        for size in [10, 50, 100]:
            frame_sampler.sample(size)
        # for each of this values, the RMSD is recalculated only for the subsample of frames and 
        # stored in the dictionary
            sampled_rmsd_vector = Property_RMSD(pro_data, frame_sampler.sampled_frame_list).rmsd
            
            sampled_rgyr_vector = Property_RadiusOfGyration(pro_data, frame_sampler.sampled_frame_list)
            
            pro_data.add_property(sampled_rmsd_vector, "RMSD", "random"+str(size))
            pro_data.add_property(sampled_rgyr_vector, "Radius Of Gyration", "random"+str(size))
            
            
            P1 = kbinDist(rmsd_vector)
            P2 = kbinDist(sampled_rmsd_vector)
            
            B_distance =   Bhatta_Distance(pro_data,P1,P2)
            
        for size in [100]:
            frame_sampler.sample(size)
        # for each of this values, the RMSD is recalculated only for the subsample of frames and 
        # stored in the dictionary
            sampled_rmsd_vector100 = Property_RMSD(pro_data, frame_sampler.sampled_frame_list).rmsd
            
            sampled_rgyr_vector100 = Property_RadiusOfGyration(pro_data, frame_sampler.sampled_frame_list)
            
            pro_data.add_property(sampled_rmsd_vector100, "RMSD", "random"+str(size))
            pro_data.add_property(sampled_rgyr_vector, "Radius Of Gyration", "random"+str(size))   
            
            
        
        
if __name__=='__main__':
    main()
