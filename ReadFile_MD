import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as LyzoData
import numpy as np
import os
def initiate_universe(topology, trajectory):
    universe = mda.Universe(topology, trajectory)
    return universe
if __name__ == '__main__':
    # File Setup
    current_directory = os.getcwd()
    
    topology = os.path.join(current_directory, '/Users/riktapatel/Desktop/Disertation/Data/MD01_1lym_example.gro')
    trajectory = os.path.join(current_directory, '/Users/riktapatel/Desktop/Disertation/Data/MD01_1lym_example.xtc')
     
LyzoData = initiate_universe(topology, trajectory)
#LyzoRef = initiate_universe(topology, trajectory)

#print(LyzoData)
type(LyzoData)
LyzoData.atoms.positions
LyzoData.trajectory
len(LyzoData.atoms)

list(LyzoData.atoms.names)
Calpha = LyzoData.select_atoms("name CA")
# len(u1.atoms)

print(Calpha)

lyzCord = np.random.rand(len(LyzoData.atoms), 3)

print(lyzCord)

# calculating RMSD between two sets of coordinates

R = rms.RMSD(LyzoData,
           select="name CA",             # superimpose on whole backbone of the whole protein
           groupselections=["name CA and (resid 1-29 or resid 60-121 or resid 160-214)",   # CORE
                            "name CA and resid 122-159",                                   # LID
                            "name CA and resid 30-59"])                                    # NMP
R.run()
R.rmsd.shape
import matplotlib.pyplot as plt
rmsd = R.rmsd.T   # transpose makes it easier for plotting
time = rmsd[1]
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)
ax.plot(time, rmsd[2], 'k-',  label="all")
ax.plot(time, rmsd[3], 'k--', label="CORE")
ax.plot(time, rmsd[4], 'r--', label="LID")
ax.plot(time, rmsd[5], 'b--', label="NMP")
ax.legend(loc="best")
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"RMSD ($\AA$)")

#for ext in ('svg', 'pdf', 'png'):
 #fig.savefig("AdK_domain_rigidity.{0}".format(ext))

fig.savefig("rmsd_5june_part2.pdf")

#ploting data 
df = pd.DataFrame(R.rmsd,
                  columns=['Frame', 'Time (ns)',
                           'name CA', 'CORE',
                           'LID', 'NMP'])

df
