import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy.linalg
import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
from sklearn import preprocessing
import os
import configparser
import random

class Protein_Data:
    
    def __init__(self,trajectory_filename,topology_filename, config_parameters):
        self.config_par = config_parameters
        self.trajectory_filename = trajectory_filename
        self.topology_filename = topology_filename
        self.trajectory_data = self._read_trajectory(self.trajectory_filename,self.topology_filename)
        self.n_frames = self.trajectory_data.trajectory.n_frames
        self.ca_atom_group = self._select_CA_atoms()
        self.property_dict = {}

    def _read_trajectory(self,trajectory_filename,topology_filename):
        trajectory_data = mda.Universe(topology_filename,trajectory_filename,permissive=False, topology_format='PDB')
        return trajectory_data

    def _select_CA_atoms(self):
        ca_atom_group = self.trajectory_data.select_atoms("name CA")
        return ca_atom_group

    def output_trj_summary(self):
        print("----------------------")
        print("TRAJECTORY INFORMATION")
        print("----------------------")
        print("n frames = {0}\nn atoms = {1}\nn CA atoms = {2}".format(self.n_frames, self.trajectory_data.trajectory.n_atoms, self.ca_atom_group.n_atoms))

    def add_property(self, property_vector, property_name, sample_label):
        if property_name in self.property_dict:
            self.property_dict[property_name][sample_label] = property_vector
        else:
            self.property_dict[property_name] = {}
            self.property_dict[property_name][sample_label] = property_vector

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
        # result.to_csv('file1.csv')
        self._statistical_analysis(result)
       
        # calling method to perfrom sampling 
        self._trajectory_sampling(result)
        # self._bhatt_distance(h)
        # Ploting graph for radius of gyration against time 
        # rgyr_df.plot(title='Radius of gyration')
        # plt.savefig(self.config_par["ImagePath"]+"rgyrgraph.png")
            
    def _statistical_analysis(self,result):
        print("----------------------")
        print("STATISTICAL ANALYSIS")
        print("----------------------")
        print(result.describe())
        # histogram for full trajectory
        result.hist()
        plt.savefig(self.config_par["ImagePath"]+"FullTraj_Hist.png")
        # scatter matrix for full trajectory
        scatter_matrix(result, alpha=0.2, figsize=(6, 6), diagonal='kde')
        plt.savefig(self.config_par["ImagePath"]+"FullTraj_Scatter.png")
    
    def _trajectory_sampling(self,result):
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
        plt.savefig(self.config_par["ImagePath"]+"Sample_Hist.png")
        # scatter matrix for subset
        scatter_matrix(subset, alpha=0.2, figsize=(6, 6), diagonal='kde')
        plt.savefig(self.config_par["ImagePath"]+"Sample_Scatter.png")

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

def main():
        config_par = get_config_parameters("SAMPLE1.INI")

        trajectory_filename = config_par["BasePath"]+config_par["trajectory"]
        topology_filename = config_par["BasePath"]+config_par["topology"]
        pro_data = Protein_Data(trajectory_filename, topology_filename, config_par)

        pro_data.output_trj_summary()

        rmsd_vector = Property_RMSD(pro_data, range(pro_data.n_frames)).rmsd
        pro_data.add_property(rmsd_vector, "RMSD", "reference")

        frame_sampler = Frame_Sampler(range(pro_data.n_frames))

        for size in [10, 50, 100]:
            frame_sampler.sample(size)

            sampled_rmsd_vector = Property_RMSD(pro_data, frame_sampler.sampled_frame_list).rmsd
            pro_data.add_property(sampled_rmsd_vector, "RMSD", "random"+str(size))

if __name__=='__main__':
    main()
