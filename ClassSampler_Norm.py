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

class Protein_Data:
    
    def __init__(self,trajectory_filename,topology_filename):
       
        self.trajectory_data = self._read_trajectory(trajectory_filename,topology_filename)    
        self.transform_data = self._transform_trajectory()
        self.normalise_data = self._normalise_trajectory()
        print("4 constructor exit")
    def _read_trajectory(self,trajectory_filename,topology_filename):
        print("2 read")
        return mda.Universe(trajectory_filename,topology_filename)
    def _transform_trajectory(self):
        print("3 transform")
    def _normalise_trajectory(self,trajectory_filename):
        
        Cov = pd.read_table( trajectory_filename, 
                    header=0, 
                    skipfooter=1,
                    delim_whitespace=True,
                    names=["A", "B", "C", "D", "E", "F", "G","H","I"]
                    )
        print(Cov)
        Dropresult = Cov.drop(['A', 'B'], axis=1)
        print(Dropresult)
        d = preprocessing.normalize(Dropresult, axis=0)
        scaled_df = pd.DataFrame(d, columns=["C", "D", "E", "F", "G","H","I"])
        scaled_df.head()
        
    def main():
        config = configparser.ConfigParser()
        config.read("SAMPLE1.INI",)
        path = config['PROTEINFILE']
        print("trajectoryfile :",path["BasePath"]+path["trajectory"])
        trajectory_filename = path["BasePath"]+path["trajectory"]
        topology_filename = path["BasePath"]+path["topology"]
   
        print(trajectory_filename)
        print(topology_filename)
        print("1 main")
        pdata1 = Protein_Data(trajectory_filename,topology_filename)
        print("object",pdata1)
        lyzo = pdata1.trajectory_data
        print("read",lyzo)
        normdata = pdata1._normalise_trajectory()
        print("Normalise",normdata)
        
    #sequence = pdata1.sequence_data
    #print(sequence)
    #lyzo_atoms = lyzo.atoms.positions
    #print("lyzo atoms",lyzo_atoms)
        calpha = lyzo.select_atoms("name CA")
        lyzo_Cord = np.random.rand(len(calpha.atoms),3)
        print(lyzo_Cord)
        lyzo.load_new(lyzo_Cord,format=MemoryReader)
        print("atom positions",lyzo.atoms.positions)
    
        print(os.getcwd())
    #glob.glob('subdir/*.pdb')

    if __name__ == '__main__':
        main()
