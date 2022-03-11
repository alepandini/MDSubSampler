"""
The utilities file consists of the following functions:

- Function that checks the memory size of the machine that the user is using.
- Function that checks the trajectory size that the user has used as an input.
- Function that checks the number of threads that are available in the system 
  (to impelent later).
- Function that checks that the XTC file is not empty and is contains at least X 
  number of frames.
- Function that checks that the number of selected residues that are used as an
  input matches the number of residues in the XTC file.
- Function that checks that the PDB/XTC files are present in the given directory.
- Funtion that checks if all given XTC files have the same number of atoms.

"""
import os
import psutil
import mdss_protein_data


def check_memory_size():
    # Getting loadover15 minutes
    # Results get updated every 5 minutes
    load1, load5, load15 = psutil.getloadavg()
    cpu_usage = (load15 / os.cpu_count()) * 100
    print("The CPU usage is: ", cpu_usage)
    print("RAM memory % used:", psutil.virtual_memory().percent)
    print(
        "The percentage of available memory is:",
        psutil.virtual_memory().available * 100 / psutil.virtual_memory().total,
    )


# Next step: Make this for a list of trajectory files
def check_trajectory_size(trajectory_file_path, topology_file_path):
    protein_data = mdss_protein_data.ProteinData(
        trajectory_file_path,
        topology_file_path,
        config_parameters=None,
    )
    trajectory_size = len(
        protein_data._frames_of_trajectory(trajectory_file_path, topology_file_path)
    )
    print("The size of the trajectory is:", trajectory_size)


# Next step: Make this for a list of trajectory files
def check_content_exists_trajectory(trajectory_file_path, topology_file_path):
    protein_data = mdss_protein_data.ProteinData(
        trajectory_file_path,
        topology_file_path,
        config_parameters=None,
    )
    if (
        len(
            protein_data._read_trajectory(
                trajectory_file_path, topology_file_path
            ).trajectory
        )
        == 0
    ):
        print("The trajectory file is empty and it has no frames")
    else:
        print(
            "The trajectory file has {} frames".format(
                len(
                    protein_data._read_trajectory(
                        trajectory_file_path, topology_file_path
                    ).trajectory
                )
            )
        )


# Testing
trajectory_file_path = "data/MD01_1lym_example_fit_short.xtc"
topology_file_path = "data/MD01_1lym_example.gro"
print("")
# print("--------------------------------")
# print("Checking the memory:")
# print("--------------------------------")
# check_memory_size()
# print("")
# print("--------------------------------")
# print("Checking the trajectory size:")
# print("--------------------------------")
# check_trajectory_size(trajectory_file_path, topology_file_path)
check_content_exists_trajectory(trajectory_file_path, topology_file_path)
