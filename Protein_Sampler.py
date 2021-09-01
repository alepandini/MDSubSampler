import MDAnalysis as mda
from MDAnalysis.analysis import rms
import MDAnalysis.analysis.pca as pca
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)
import numpy.linalg
import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix
from matplotlib import pyplot
from sklearn import preprocessing
import os
import configparser
import random
import dictances 
from dictances import bhattacharyya, bhattacharyya_coefficient,kullback_leibler
from scipy.stats import ks_2samp
import itertools 
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu
import pingouin as pg
import csv


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
        
    
    def statistical_analysis(self,stat_vector,stat_sampled_vector,sample_label):
        
        print("----------------------")
        print("STATISTICAL ANALYSIS")
        print("----------------------")
        vector_mean = np.mean(np.array([stat_vector]), axis=0)
        vector_variance = np.var(stat_vector)
        
        print("mean for {0} vector::{1}".format(sample_label,vector_mean))
        print("variance for {0} vector::{1}".format(sample_label,vector_variance)) 
        
        print("----------------------")
        print("Kolmogorov-Smirnov test")
        print("----------------------")
       
        print(ks_2samp(stat_vector, stat_sampled_vector))
        
        print(stats.mannwhitneyu(stat_vector, stat_sampled_vector))
        print(pg.ttest(x=stat_vector,y=vector_mean))
        print(pg.ttest(x=stat_vector, y=stat_sampled_vector, correction=False).round(2))
        print(pg.ttest(x=stat_sampled_vector, y=stat_vector, correction=False).round(2))

        
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
       
        self.rgyr = []
        self.time = []
        for frame in frame_list:
            protein_data.trajectory_data.trajectory[frame]
            self.time.append(protein_data.trajectory_data.trajectory.time)
            self.rgyr.append(protein_data.trajectory_data.select_atoms(atom_selection).radius_of_gyration())
        # ploting radius of Gyration and saving as PDF in root directory
        #ax = plt.subplot(111)
        #ax.plot(self.time,self.rgyr, 'b--', lw=2, label=r"$R_G$")
        #ax.set_xlabel("time (ps)")
        #ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)")
        #ax.figure.savefig("Rgyr.pdf")
        #plt.draw()

class Bhatta_Distance:
       
        def __init__(self, Prop_vector, verbose=True):
            
            min_value = np.min(Prop_vector)
            max_value = np.max(Prop_vector)

            # discretisation of the original vector with all values
            freq_prop_vector = discretize_to_dict(Prop_vector, min_value, max_value)
         
            # test - the distance should be zero (or close to zero) on itself
            sample_size = len(Prop_vector)
            b_distance = dictances.bhattacharyya(freq_prop_vector, freq_prop_vector)
            print("size: {0:4d} distance: {1:.2f}".format(sample_size, b_distance))

            print("-------------------------")         

            for sample_size in [10,20,50,100,200,500]:
            
                sub_prop_vector = random.sample(Prop_vector, sample_size)

            # discretisation of the subsampled vector
                freq_sub_prop_vector = discretize_to_dict(sub_prop_vector, min_value, max_value)
                
                b_distance = dictances.bhattacharyya(freq_prop_vector, freq_sub_prop_vector)
                print("size: {0:4d} distance: {1:.2f}".format(sample_size, b_distance))
                
            
class KL_diver:
        def __init__(self, Prop_vector, verbose=True):
            
            min_value = np.min(Prop_vector)
            max_value = np.max(Prop_vector)

            # discretisation of the original vector with all values
            KL_freq_prop_vector = discretize_to_dict(Prop_vector, min_value, max_value)
            
            #  replace 0 with small number and then rescale values to make the sum equal to 1           
            KL_freq_prop_vector_clean = replace_zero(KL_freq_prop_vector)
            
            # test - the distance should be zero (or close to zero) on itself
            sample_size = len(Prop_vector)
            
            kl_pq_distance = dictances.kullback_leibler(KL_freq_prop_vector_clean, KL_freq_prop_vector_clean)
            
            print("size: {0:4d} distance: {1:.2f}".format(sample_size, kl_pq_distance))
            
            print("-------------------------")         

            for sample_size in [10,20,50,100,200,500]:
            
                KL_sub_prop_vector = random.sample(Prop_vector, sample_size)

                # discretisation of the subsampled vector
                KL_freq_sub_prop_vector = discretize_to_dict(KL_sub_prop_vector, min_value, max_value)
                
                KL_freq_sub_prop_vector_clean = replace_zero(KL_freq_sub_prop_vector)
                # calculate the kl divergence
                kl_pq_distance = dictances.kullback_leibler(KL_freq_sub_prop_vector_clean, KL_freq_sub_prop_vector_clean)
                
                print("size: {0:4d} distance: {1:.2f}".format(sample_size, kl_pq_distance))
          
        
class PCA_analysis:
    
        def __init__(self,protein_data, frame_list, atom_selection = "name CA", verbose=True):
           
               
            PSF_pca = pca.PCA(protein_data.trajectory_data, select='name CA',align=False, mean=None).run()
               
            n_pcs = np.where(PSF_pca.results.cumulated_variance > 0.95)[0][0]
               
            pca_space = PSF_pca.transform(protein_data.trajectory_data.select_atoms('name CA'), n_components=n_pcs)
                
            print("----------------------")
            print("PCA")
            print("----------------------")
            print("PCA Shape:",pca_space.shape)
            print("PCA Variance:",PSF_pca.results.cumulated_variance[0])
            print("PCA Variance:",PSF_pca.results.cumulated_variance[2])
            #   reduce component to 5              
                
            transformed = PSF_pca.transform(protein_data.trajectory_data.select_atoms('name CA'),
                                                n_components=5)
                
            print("Reduced PCA component shape::",transformed.shape)
                
            df = pd.DataFrame(transformed,columns=['PC{}'.format(i+1) for i in range(5)])
            df['Time (ps)'] = df.index * protein_data.trajectory_data.trajectory.dt
            print(df.head())
                
                
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

def discretize_to_dict(values, min_value, max_value):
        bin_size = (max_value - min_value) / 100.0
        bin_vector = np.arange(min_value, max_value, bin_size)
        counts, bins  = np.histogram(values, bins = bin_vector)
        freq_dict = dict(zip(bins, counts/len(values)))
        return(freq_dict)
def calculate_statistic(rmsd_vector):
        summary_stat = np.array(rmsd_vector)
        res = np.mean(summary_stat)
        return summary_stat

def replace_zero(freq_vector):
        
        dic_out = {}
        for x, y in freq_vector.items():
            if y != 0:
                dic_out[x] = y 
            else:
                dic_out[x] = 0.0001
    
        value_sum = sum(dic_out.values())
        dic_out = {k: v / value_sum for k, v in dic_out.items()}
        return dic_out 
    
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
        
        
        print("-----BHATTACHARYA DISTANCE FOR RMSD--------------------")  
        b_dict = Bhatta_Distance(rmsd_vector)
        print("-----BHATTACHARYA DISTANCE FOR RGYR--------------------")  
        b_rgvr_dict = Bhatta_Distance(rgyr_vector.rgyr)
        print("-----KL divergence DISTANCE FOR RMSD--------------------") 
        KL_Dict_measure = KL_diver(rmsd_vector)
        print("-----KL divergence DISTANCE FOR RGYR--------------------")  
        b_rgvr_dict = KL_diver(rgyr_vector.rgyr)
        
        
        pro_data.add_property(rmsd_vector, "RMSD", "reference")
        pro_data.add_property(rgyr_vector, "Radius Of Gyration", "reference")
        #       creates a Frame_Sampler object
        frame_sampler = Frame_Sampler(range(pro_data.n_frames))
        #       for different values of sample size, the sampler randomly selects frames

        for size in [10,20,50,100,200,500]:
            frame_sampler.sample(size)
            
        # for each of this values, the RMSD is recalculated only for the subsample of frames and 
        # stored in the dictionary
            sampled_rmsd_vector = Property_RMSD(pro_data, frame_sampler.sampled_frame_list).rmsd
            
            sampled_rgyr_vector = Property_RadiusOfGyration(pro_data, frame_sampler.sampled_frame_list)
            
            pro_data.add_property(sampled_rmsd_vector, "SAMPLED RMSD", "random"+str(size))
            pro_data.add_property(sampled_rgyr_vector, "SAMPLED Radius Of Gyration", "random"+str(size))
            
        pro_data.statistical_analysis(rmsd_vector,sampled_rmsd_vector,"RMSD")    
        pro_data.statistical_analysis(rgyr_vector.rgyr,sampled_rgyr_vector.rgyr,"RGYR")
       
        
        pca_vector = PCA_analysis(pro_data, range(pro_data.n_frames))
if __name__=='__main__':
    main()
