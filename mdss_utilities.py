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
- Function that checks if all given XTC files have the same number of atoms.

"""
import os
import psutil
import mdss_protein_data


def check_memory_size(trajectory_file_path, topology_file_path):
    mem = psutil.virtual_memory()
    file_size_trajectory = os.path.getsize(trajectory_file_path)
    file_size_topology = os.path.getsize(topology_file_path)
    if file_size_trajectory > mem.available:
        raise Exception(
            "The size of the trajectory file is larger than the available memory"
        )
    if file_size_topology > mem.available:
        raise Exception(
            "The size of the topology file is larger than the available memory"
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
    print("The trajectory has {} number of frames:".format(trajectory_size))


# for loop
# Next step: Make this for a list of trajectory files
def check_content_exists_trajectory(trajectory_file_path, topology_file_path):
    protein_data = mdss_protein_data.ProteinData(
        trajectory_file_path,
        topology_file_path,
        config_parameters=None,
    )
    # initialise the len
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


# add error (exception file not exist)
def check_files_exist(trajectory_file_path, topology_file_path):
    if not os.path.isfile(trajectory_file_path):
        print("The XTC trajectory file was not found in the directory")
    if not os.path.isfile(topology_file_path):
        print("The PDB topology file was not found in the directory")


# This will run only if this file is run as a script
if __name__ == "__main__":
    trajectory_file_path = "data/MD01_1lym_example_fit_short.xtc"
    topology_file_path = "data/MD01_1lym_example.gro"
    check_memory_size(trajectory_file_path, topology_file_path)
    # check_trajectory_size(trajectory_file_path, topology_file_path)
    # check_content_exists_trajectory(trajectory_file_path, topology_file_path)
    # check_files_exist(trajectory_file_path, topology_file_path)
