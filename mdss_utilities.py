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


class NotEnoughMemoryError(Exception):
    pass


class EmptyTrajectoryError(Exception):
    pass


def check_file_size(filepath):
    mem = psutil.virtual_memory()
    file_size = os.path.getsize(filepath)
    if file_size > mem:
        raise NotEnoughMemoryError(
            "The size of {} is larger than the available memory".format(filepath)
        )


def check_trajectory_size(trajectory_file_path, topology_file_path):
    protein_data = mdss_protein_data.ProteinData(
        trajectory_file_path,
        topology_file_path,
        config_parameters=None,
    )
    trajectory_size = len(
        protein_data._read_trajectory(
            trajectory_file_path, topology_file_path
        ).trajectory
    )
    if trajectory_size == 0:
        raise EmptyTrajectoryError(
            "The trajectory has no frames in it "
            "(trajectory file = {}, "
            "topology file = {})".format(trajectory_file_path, topology_file_path)
        )


def check_multiple_trajectories_size(list_of_traj, list_of_top):
    errors = []
    for traj, top in zip(list_of_traj, list_of_top):
        try:
            check_trajectory_size(traj, top)
        except EmptyTrajectoryError as e:
            errors.add(e)

    if not errors:
        return

    for error in errors:
        print(error)
        raise RuntimeError()


def check_file_exists(filepath):
    if not os.path.isfile(filepath):
        raise FileNotFoundError("File {} does not exist".format(filepath))


# This will run only if this file is run as a script
if __name__ == "__main__":
    trajectory_file_path = "data/MD01_1lym_example_fit_short.xtc"
    topology_file_path = "data/MD01_1lym_example.gro"
    list_of_trajectories = [
        "data/MD01_1lym_example_fit_short.xtc",
        "data/MD01_1lym_example_fit_short.xtc",
    ]
    list_of_topologies = ["data/MD01_1lym_example.gro", "data/MD01_1lym_example.gro"]
    check_file_size(trajectory_file_path, topology_file_path)
    check_trajectory_size(trajectory_file_path, topology_file_path)
    check_file_exists(trajectory_file_path, topology_file_path)
    check_multiple_trajectories_size(list_of_trajectories, list_of_topologies)
