'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Rikta Patel

    
     This file is part of the prototype framework development of subsampler for large biomolecular trajectories.  
     (https://github.com/alepandini/MDSubSampler).
     Copyright (c) 2020-21 Rikta Patel and Alessandro Pandini.

     This program is free software: you can redistribute it and/or modify 
     it under the terms of the GNU General Public License as published by  
     the Free Software Foundation, version 3.

     This program is distributed in the hope that it will be useful, but 
     WITHOUT ANY WARRANTY; without even the implied warranty of 
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
     General Public License for more details.

     You should have received a copy of the GNU General Public License 
     along with this program. If not, see <http://www.gnu.org/licensÅes/>.

'''
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import MDAnalysis.analysis.pca as pca
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)
import numpy.linalg
import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix
from matplotlib import pyplot
import matplotlib.pyplot as plt
from sklearn import preprocessing
import os
import configparser
import random
import dictances 
from dictances import bhattacharyya, bhattacharyya_coefficient,kullback_leibler
import itertools 
import seaborn as sns
import scipy
from scipy import stats
from scipy.stats import mannwhitneyu, shapiro, ks_2samp, ttest_ind
from scipy.signal import savgol_filter
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
        trajectory_data = mda.Universe(topology_filename,trajectory_filename,permissive=False,
                                       topology_format='PDB')
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
        print("-----------MAN-WI-----------")
        statMW, pMW = mannwhitneyu(stat_vector, stat_sampled_vector)
        print('stat=%.3f, p=%.3f' % (statMW, pMW))
        if pMW > 0.05:
            print('Probably the same distribution')
        else:
            print('Probably different distributions')
#         print(pg.ttest(x=stat_vector, y=stat_sampled_vector, correction=True).round(2))
        print("----------T-TEST ONE------------")
        stat1, p1 = ttest_ind(stat_vector, stat_sampled_vector)
        print('stat1=%.3f, p1=%.3f' % (stat1, p1))
        if p1 > 0.05:
            print('Probably the same distribution')
        else:
            print('Probably different distributions')
      
        print("----------------------")
        print("Shapio test")
        print("----------------------")
        stat, p = shapiro(stat_vector)
        print('stat=%.3f, p=%.3f' % (stat, p))
        if p > 0.05:
            print('Probably Gaussian')
        else:
            print('Probably not Gaussian')
        
        
    def add_property(self, property_vector, property_name, sample_label):
        if property_name in self.property_dict:
            self.property_dict[property_name][sample_label] = property_vector
            
        else:
            self.property_dict[property_name] = {}
            self.property_dict[property_name][sample_label] = property_vector

class Frame_Sampler:
        def __init__(self, frame_list, seed_number = 1999):
            random.seed(seed_number)
            self.frame_list = frame_list    
            self.sampled_frame_list = None

        def sample(self, size):
            self.sampled_frame_list = random.sample(self.frame_list, size)
            self.sampled_frame_list.sort()
#         def stat_sample(self,size):
            
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
        
class Bhatta_Distance:
       
        def __init__(self, Prop_vector, verbose=True):
            
            min_value = np.min(Prop_vector)
            max_value = np.max(Prop_vector)
            print("B  size",len(Prop_vector))
            # discretisation of the original vector with all values
            freq_prop_vector = discretize_to_dict(Prop_vector, min_value, max_value)
            print("B disatacne size",len(freq_prop_vector))
            # test - the distance should be zero (or close to zero) on itself
            sample_size = len(freq_prop_vector)
            b_distance = dictances.bhattacharyya(freq_prop_vector, freq_prop_vector)
            print("size: {0:6d} distance: {1:.3f}".format(sample_size, b_distance))

            print("-------------------------")         

            for sample_size in [10,50,100,200,400,500]:
            #for sample_size in [800,4000,8000,16000,32000,40000]:
            
                sub_prop_vector = random.sample(Prop_vector, sample_size)

            # discretisation of the subsampled vector
                freq_sub_prop_vector = discretize_to_dict(sub_prop_vector, min_value, max_value)
                
                b_distance = dictances.bhattacharyya(freq_prop_vector, freq_sub_prop_vector)
                print("size: {0:6d} distance: {1:.3f}".format(sample_size, b_distance))
                
            
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
            
            print("size: {0:6d} distance: {1:.3f}".format(sample_size, kl_pq_distance))
            
            print("-------------------------")         

            for sample_size in [10,50,100,200,400,500]:
            #for sample_size in [800,4000,8000,16000,32000,40000]:
            
                KL_sub_prop_vector = random.sample(Prop_vector, sample_size)
                
                # discretisation of the subsampled vector
                KL_freq_sub_prop_vector = discretize_to_dict(KL_sub_prop_vector, min_value, max_value)
                
                KL_freq_sub_prop_vector_clean = replace_zero(KL_freq_sub_prop_vector)
                # calculate the kl divergence
                kl_pq_distance = dictances.kullback_leibler(KL_freq_sub_prop_vector_clean, KL_freq_sub_prop_vector_clean)
                
                print("size: {0:6d} distance: {1:.3f}".format(sample_size, kl_pq_distance))
class Pearson_dict:
        def __init__(self, Prop_vector, verbose=True):
            
            min_value = np.min(Prop_vector)
            max_value = np.max(Prop_vector)

            # discretisation of the original vector with all values
            KL_freq_prop_vector = discretize_to_dict(Prop_vector, min_value, max_value)
            
            #  replace 0 with small number and then rescale values to make the sum equal to 1           
            KL_freq_prop_vector_clean = replace_zero(KL_freq_prop_vector)
            
            # test - the distance should be zero (or close to zero) on itself
            sample_size = len(Prop_vector)
            
            kl_pq_distance = dictances.pearson(KL_freq_prop_vector, KL_freq_prop_vector)
            
            print("size: {0:6d} distance: {1:.3f}".format(sample_size, kl_pq_distance))
            
            print("-------------------------")         

            for sample_size in [10,50,100,200,400,500]:
            #for sample_size in [800,4000,8000,16000,32000,40000]:
            
                KL_sub_prop_vector = random.sample(Prop_vector, sample_size)

                # discretisation of the subsampled vector
                KL_freq_sub_prop_vector = discretize_to_dict(KL_sub_prop_vector, min_value, max_value)
                
                KL_freq_sub_prop_vector_clean = replace_zero(KL_freq_sub_prop_vector)
                # calculate the kl divergence
                kl_pq_distance = dictances.pearson(KL_freq_sub_prop_vector, KL_freq_sub_prop_vector)
                
                print("size: {0:6d} distance: {1:.3f}".format(sample_size, kl_pq_distance))          
        

class Property_PCA_analysis:
    
        def __init__(self,protein_data, frame_list, atom_selection = "name CA", verbose=True):
           
            calpha_pca = protein_data.trajectory_data.select_atoms('name CA')    
            full_pca = pca.PCA(protein_data.trajectory_data, select='name CA',align=False, mean=None).run()
            print("PCA TYPE ", type(full_pca))
            full_trans = full_pca.transform(calpha_pca,n_components=5)
            df1 = pd.DataFrame(full_trans ,columns=['PC{}'.format(i+1) for i in range(5)])
            df1['Time (ps)'] = df1.index * protein_data.trajectory_data.trajectory.dt
            print(df1.head())
            print("Variance",full_pca.cumulated_variance[5])
            plt.plot(full_pca.cumulated_variance[:5])
           
            # pcgenerate your frames array by any means necessary, this is an example
            
class Subsampled_PCA_analysis:
     def __init__(self,protein_data, frame_list, atom_selection = "name CA", verbose=True):
            calpha_pca = protein_data.trajectory_data.select_atoms('name CA')
            full_pca = pca.PCA(protein_data.trajectory_data, select='name CA',align=False, mean=None).run()
#             frames = np.array([10,50,100,200,400,500,1000])   
            sliced_traj = protein_data.trajectory_data.trajectory[frame_list]
            print("Slice array",sliced_traj)
            coordinates = np.empty((len(sliced_traj), protein_data.trajectory_data.select_atoms('name CA').n_atoms, 3), dtype=np.float32)
            for i, ts in enumerate(sliced_traj):
                coordinates[i] = protein_data.trajectory_data.select_atoms(atom_selection).positions
            u2 = mda.Merge(calpha_pca)            # create the caplha-only Universe
            u2.load_new(coordinates, format=MemoryReader)
            # the u2 universe now contains the c-alpha with only the frames of interest
            # use the subsampled universe u2 (should also be faster because its in memory)
            ca = u2.select_atoms("name CA")
            sampled_pca = pca.PCA(u2, select='name CA',align=False, mean=None).run()
            print("PCA SAMPLED TYPE ", type(sampled_pca))
            sampled_transformed = sampled_pca.transform(ca,n_components=5)
            
            df = pd.DataFrame(sampled_transformed ,columns=['PC{}'.format(i+1) for i in range(5)])
            df['Time (ps)'] = df.index * protein_data.trajectory_data.trajectory.dt
            print(df.head())
#             print("Variance",sampled_transformed.variance[0])
            print("Variance subsampled ",sampled_pca.cumulated_variance[5])
            print("PCA GRAPH") 
            sns.lmplot( x="PC1", y="PC2",
                           data=df, 
                           fit_reg=False, 
                          legend=True
                         ) 
            plt.show()

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
def Property_RMSF(protein_data, frame_list):
        protein_data.trajectory_data.trajectory[0]
    
        rmsf = []
        for frame in frame_list:
            protein_data.trajectory_data.trajectory[frame]
            R = rms.RMSF(protein_data.trajectory_data.select_atoms('name CA')).run()
            return R.rmsf

def plot_histogram(hist_vector):
            
            noise_filter = savgol_filter(hist_vector, 5, 2) 
            _, bins, _ = plt.hist(noise_filter, 20, density=1, alpha=0.5)
           
            mu, sigma = scipy.stats.norm.fit(noise_filter)
            best_fit_line = scipy.stats.norm.pdf(bins, mu, sigma)
            plt.ylim([0, 1.4])
            plt.xlim([0.5, 4.0])
            plt.plot(bins, best_fit_line)
            plt.show()
def plot_boxplot(box_vector):
            print("Box plot")
            plt.boxplot(box_vector)
            plt.show()
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
    
def save_output_plot(RGYR_result):
    stat_method_dict = {
          'RGYR': RGYR_result,
          #'RMSD_CAlpha': RMSD_CAlpha_statistic,
          #'RMSD_Full_Protein': RMSD_Full_Protein_statistic,
        #'RMSF': RMSF_statistic,
            }
    method = "RGYR"
    i = [10,20,50,100,200,500]
    for (method, statistic) in stat_method_dict.items():
      
        for size in i:
            ax = plt.subplot(111)
            W = savgol_filter(RGYR_result.time, 5, 2) 
            W1 = savgol_filter(RGYR_result.rgyr, 5, 2) 
            ax.plot(W,W1, 'steelblue', lw=2, label=r"$R_G$")
            ax.set_xlabel("time (ps)")
            ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)")
            ax.set_ylim([15.75, 17.50])
        # ax.set_xlim([500, 2500])
        #ax.figure.savefig("Rgyrsample.pdf")
        plt.show()
        plt.savefig(method+str(size)+'.png')
       
    
def main():
        config_par = get_config_parameters("SAMPLE1.INI")
        
        trajectory_filename = config_par["BasePath"]+config_par["trajectory"]
        topology_filename = config_par["BasePath"]+config_par["topology"]
        #       creates a Protein_Data object
        pro_data = Protein_Data(trajectory_filename, topology_filename, config_par)
       
        pro_data.output_trj_summary()
        
        #       calculates an example property
        rmsd_vector = Property_RMSD(pro_data, range(pro_data.n_frames)).rmsd
        
        hist_rmsd = plot_histogram(rmsd_vector)
         #box_rmsd = plot_boxplot(rmsd_vector)
         
        rmsf_vector = Property_RMSF(pro_data, range(pro_data.n_frames))
        
        rgyr_vector = Property_RadiusOfGyration(pro_data, range(pro_data.n_frames))
        
        
        #ploting radius of Gyration and saving as PDF in root directory
        
        print("-----BHATTACHARYA DISTANCE FOR RMSD--------------------")  
        b_dict = Bhatta_Distance(rmsd_vector)
#         print("-----BHATTACHARYA DISTANCE FOR RMSF--------------------")  
#         b_rmsf_dict = Bhatta_Distance(rmsf_vector)
        print("-----BHATTACHARYA DISTANCE FOR RGYR--------------------")  
        b_rgvr_dict = Bhatta_Distance(rgyr_vector.rgyr)
        print("-----Kullback–Leibler divergence FOR RMSD--------------------") 
        KL_Dict_measure = KL_diver(rmsd_vector)
        print("-----Kullback–Leibler divergence DISTANCE FOR RGYR--------------------")  
        b_rgvr_dict = KL_diver(rgyr_vector.rgyr)
        print("-----PEARSON DISTANCE FOR RMSD--------------------") 
        pearson_dict = Pearson_dict(rmsd_vector)
        print("-----PEARSON DISTANCE DISTANCE FOR RGYR--------------------")  
        pearson_rgvr_dict = Pearson_dict(rgyr_vector.rgyr)
        print("-----PCA--------------------") 
        
        pca_vector = Property_PCA_analysis(pro_data, range(pro_data.n_frames))
         #pca_vector1 = Property_PCA_analysis1(pro_data, range(pro_data.n_frames))
        
        pro_data.add_property(rmsd_vector, "RMSD", "reference")
        # pro_data.add_property(rmsf_vector, "RMSF", "reference")
        pro_data.add_property(rgyr_vector, "Radius Of Gyration", "reference")
        pro_data.add_property(pca_vector, "PCA", "reference")
        #       creates a Frame_Sampler object
        frame_sampler = Frame_Sampler(range(pro_data.n_frames))
        #       for different values of sample size, the sampler randomly selects frames

        for size in [10,50,100,200,400,500]:
        
        #for size in [800,4000,8000,16000,32000,40000]:
            frame_sampler.sample(size)
            
        # for each of this values, the RMSD is recalculated only for the subsample of frames and 
        # stored in the dictionary
            sampled_rmsd_vector = Property_RMSD(pro_data, frame_sampler.sampled_frame_list).rmsd
        
            
            hist_sampled_rmsd = plot_histogram(sampled_rmsd_vector)
            #box_sampled_rmsd = plot_boxplot(sampled_rmsd_vector)
            
            #sampled_rmsf_vector = Property_RMSF_sampled(pro_data, frame_sampler.sampled_frame_list)
           
            sampled_rgyr_vector = Property_RadiusOfGyration(pro_data,frame_sampler.sampled_frame_list)
             
            save_output_plot(sampled_rgyr_vector)
            
            sampled_pca_vector = Subsampled_PCA_analysis(pro_data, frame_sampler.sampled_frame_list)
            
            pro_data.add_property(sampled_rmsd_vector, "SAMPLED RMSD", "random"+str(size))
            #pro_data.add_property(sampled_rmsf_vector, "SAMPLED RMSF", "random"+str(size))
            pro_data.add_property(sampled_rgyr_vector, "SAMPLED Radius Of Gyration", "random"+str(size))
#             pro_data.add_property(sampled_pca_vector, "SAMPLED PCA", "random"+str(size))
            
        pro_data.statistical_analysis(rmsd_vector,sampled_rmsd_vector,"RMSD")    
        pro_data.statistical_analysis(rgyr_vector.rgyr,sampled_rgyr_vector.rgyr,"RGYR")
        
if __name__=='__main__':
    main()
