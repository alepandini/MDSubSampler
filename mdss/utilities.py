"""
The utilities file consists of the following functions:

- Function that checks the memory size of the machine that the user is using.
- Function that checks the trajectory size that the user has used as an input.
- Function that checks the number of threads that are available in the system 
  (to impelent later).
- Function that checks that the XTC and topology files are not empty and contain 
  at least X number of frames.
- Function that checks that the number of selected residues that are used as an
  input matches the number of residues in the XTC file.
- Function that checks that the PDB/XTC files are present in the given directory.

"""
import os
import psutil
import mdss.protein_data as pd
import mdss.graph as g


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
    protein_data = pd.ProteinData(
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
    return trajectory_size


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


def check_number_of_residues(trajectory_file_path, topology_file_path, atom_selection):
    protein_data = pd.ProteinData(
        trajectory_file_path,
        topology_file_path,
        config_parameters=None,
    )
    if len(atom_selection) > len(
        protein_data.trajectory_data.select_atoms(atom_selection)
    ):
        raise Exception("The atom selection is greater than the imported XTC file")


def check_file_exists(filepath):
    if not os.path.isfile(filepath):
        raise FileNotFoundError("File {} does not exist".format(filepath))


def write_output_files(output_folder, file_prefix, p_prop, s_prop, p_data, p=None):
    p_format = "_" if p is None else f"_{int(p)}_"

    filename = "{}{}{}.dat".format(
        file_prefix,
        p_format,
        p_prop.display_name,
    )
    filepath = os.path.join(output_folder, filename)

    s_prop.write_property_vector(filepath)

    filename = "{}{}{}.xtc".format(
        file_prefix,
        p_format,
        p_prop.display_name,
    )
    filepath = os.path.join(output_folder, filename)
    selected_frames = p_data.frame_selection_indices(s_prop.frame_indices)
    p_data.write_xtc_file(filepath, selected_frames)

    filename = "{}{}{}".format(
        file_prefix,
        p_format,
        p_prop.display_name,
    )
    filepath = os.path.join(output_folder, filename)
    subsampled_traj = p_data.frame_selection_iterator(s_prop.frame_indices)
    p_data.cast_output_traj_to_numpy(filepath, subsampled_traj)


def plot_property(output_folder, file_prefix, p_prop, s_prop, p=None):
    p_format = "_" if p is None else f"_{int(p)}_"
    filename = "{}{}{}_{}.png".format(
        file_prefix, p_format, p_prop.display_name, "plot"
    )
    filepath = os.path.join(output_folder, filename)
    graph = g.PropertyPlot(p_prop, s_prop, filepath)
    graph.plot(p_prop.display_name, p_prop, s_prop, filepath)


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
    check_number_of_residues(trajectory_file_path, topology_file_path, "name CA")
