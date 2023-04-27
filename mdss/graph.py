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
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math


class PropertyPlot:
    """
    Represents vector plot after property calculation

    Attributes
    -----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    outfilepath : str
        path where output file with plot will be saved
    """

    def __init__(self, property, sampled_property, outfilepath):
        self.property = property
        self.property_vector = property.property_vector
        self.sampled_property = sampled_property
        self.sampled_property_vector = sampled_property.property_vector
        self.frame_indices = property.frame_indices
        self.outfilepath = outfilepath

    def convert_list_to_data_frame(self, calculated_property):
        """
        Converts list with calculated property to pandas dataframe

        Attributes
        -----------
        calculated_property: list
            A list that contains a specific calculated property for protein trajectory
        """
        df = pd.DataFrame(calculated_property)
        return df

    def plot(
        self, property_name, reference_df, sample_df, sample_size, outfilepath=None
    ):
        """
        Plots overlapped distribution of full and sample trajectory and saves file

        Attributes
        -----------
        property_name: str
            Name of property that is calculated
        reference_df: dataframe
            dataframe with calculated property for full protein trajectory
        sample_df: dataframe
            dataframe with calculated property for sample protein trajectory
        sample_size: int
            size of sample trajectory
        outfilepath : str
            Path where output file with plot will be saved
        """
        if outfilepath is None:
            outfilepath = self.outfilepath

        reference_df = self.convert_list_to_data_frame(self.property.property_vector)
        sample_df = self.convert_list_to_data_frame(
            self.sampled_property.property_vector
        )

        min_plot_value = min(float(reference_df.min()), float(sample_df.min()))
        max_plot_value = max(float(reference_df.max()), float(sample_df.max()))
        n_breaks = 50

        breaks_seq = []
        for i in np.arange(
            min_plot_value, max_plot_value, (max_plot_value - min_plot_value) / n_breaks
        ):
            breaks_seq.append(i)

        plt.style.use("ggplot")
        plt.hist(reference_df, breaks_seq, alpha=0.5, label="reference", density=True)
        plt.hist(sample_df, breaks_seq, alpha=0.5, label="sample", density=True)
        plt.xlim(math.floor(min_plot_value), math.ceil(max_plot_value))
        plt.ylim(0, 1)
        plt.xlabel(
            "{}{}".format(property_name, "/Å"),
            fontsize="12",
            horizontalalignment="center",
        )
        plt.ylabel(
            "density",
            fontsize="12",
            horizontalalignment="center",
        )
        plt.legend(loc="upper right")
        plt.title(
            "{} {} {} {} {} {}".format(
                "Sample size:",
                sample_size,
                "\n",
                self.sampled_property.dissimilarity_name,
                "distance:",
                round(self.sampled_property.ref_dissimilarity, 5),
            )
        )

        plt.savefig(outfilepath)
        plt.clf()
