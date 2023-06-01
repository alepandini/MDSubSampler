"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Namir Oues
    
    This file is part of the MDSubSampler software 
    (https://github.com/alepandini/MDSubSampler).
    Copyright (c) 2023 Namir Oues and Alessandro Pandini.

    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
import os
import psutil
import mdss.graph as g
from mdss.log_setup import log
import mdss.protein_data as pd


class NotEnoughMemoryError(Exception):
    pass


class EmptyTrajectoryError(Exception):
    pass


def check_file_size(filepath):
    """
    Checks file size and machine memory to ensure there is enough memory in the system
    """
    mem = psutil.virtual_memory()
    file_size = os.path.getsize(filepath)
    if file_size > mem:
        log.error(
            "{:15s} The size of {} is larger than the available memory".format(
                "STEPS", filepath
            )
        )
        raise NotEnoughMemoryError(
            "The size of {} is larger than the available memory".format(filepath)
        )


def check_trajectory_size(trajectory_file_path, topology_file_path):
    """
    Checks trajectory size and that trajectory and topology files are not empty and
    contain at least X number of frames.
    """
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
        log.error("{:15s} The trajectory has no frames in it ".format("STEPS"))
        raise EmptyTrajectoryError(
            "The trajectory has no frames in it "
            "(trajectory file = {}, "
            "topology file = {})".format(trajectory_file_path, topology_file_path)
        )
    return trajectory_size


def check_multiple_trajectories_size(list_of_trajectory_files, topology_file_path):
    """
    Checks size of multiple trajectories and that trajectory and topology files are not
    empty and contain at least X number of frames.
    """
    errors = []
    for traj in list_of_trajectory_files:
        try:
            check_trajectory_size(traj, topology_file_path)
        except EmptyTrajectoryError as e:
            errors.add(e)

    if not errors:
        return

    for error in errors:
        print(error)
        raise RuntimeError()


def check_number_of_residues(trajectory_file_path, topology_file_path, atom_selection):
    """
    Checks that the number of selected residues that are used as an
    input matches the number of residues in the XTC file.
    """
    protein_data = pd.ProteinData(
        trajectory_file_path,
        topology_file_path,
        config_parameters=None,
    )
    if len(atom_selection) > len(
        protein_data.trajectory_data.select_atoms(atom_selection)
    ):
        log.error(
            "{:15s} The atom selection is greater than the imported XTC file".format(
                "STEPS"
            )
        )
        raise Exception("The atom selection is greater than the imported XTC file")


def check_file_exists(filepath):
    """
    Checks that file exists in a given filepath
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError("File {} does not exist".format(filepath))


def write_output_files(
    output_folder,
    file_prefix,
    p_prop,
    s_prop,
    p_data,
    p=None,
    unit="nanometer",
    machine_learning=False,
):
    """
    Writes all output files

    Attributes
    -----------
    output_folder: str,
            file path for output folder given as user input
    file_prefix: str,
            prefix given as user input
    p_prop: list
            list with values of property for full trajectory to plot
    s_prop: list
            list with values of property for sample trajectory to plot
    p_data: ProteinData Class object

    p: int
       percentage of sample trajectory
    unit: str
       unit that the property will be calculated and saved
    """
    if output_folder is None:
        return

    p_format = "_" if p is None else f"_{p}_"

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
    coordinates_array = p_data.cast_output_traj_to_numpy(
        filepath, subsampled_traj, unit
    )
    if machine_learning:
        filename = "{}{}{}{}".format(
            file_prefix, p_format, p_prop.display_name, "_ML_input"
        )
        filepath = os.path.join(output_folder, filename)
        ML_input = p_data.convert_numpy_to_2D(coordinates_array, filepath)
        filename_train = "{}{}{}{}".format(
            file_prefix, p_format, p_prop.display_name, "_ML_train"
        )
        filepath_train = os.path.join(output_folder, filename_train)
        filename_test = "{}{}{}{}".format(
            file_prefix, p_format, p_prop.display_name, "_ML_test"
        )
        filepath_test = os.path.join(output_folder, filename_test)
        p_data.ML_input_prep(ML_input, filepath_train, filepath_test)


def plot_property(output_folder, file_prefix, p_prop, s_prop, p=None):
    """
    Plots overlapped property distributions of full and sample trajectory

    Attributes
    -----------
    output_folder: str,
            file path for output folder given as user input
    file_prefix: str,
            prefix that was given as a choice by the user
    p_prop: list
            list with values of property for full trajectory to plot
    s_prop: list
            list with values of property for sample trajectory to plot
    p: int
       percentage of sample trajectory

    Returns
    -----------
    png file with overlay property distribution of full and sample trajectories

    """
    if output_folder is None:
        return

    p_format = "_" if p is None else f"_{p}_"
    filename = "{}{}{}_{}.png".format(
        file_prefix, p_format, p_prop.display_name, "plot"
    )
    filepath = os.path.join(output_folder, filename)
    graph = g.PropertyPlot(p_prop, s_prop, filepath)
    return graph.plot(p_prop.display_name, p_prop, s_prop, p, filepath)
