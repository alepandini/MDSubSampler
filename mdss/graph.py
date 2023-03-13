import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math


class PropertyPlot:
    """
    Represents the vector plot after property calculation

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
        self.property_name = property.display_name
        self.sampled_property_vector = sampled_property.property_vector
        self.frame_indices = property.frame_indices
        self.outfilepath = outfilepath

    def convert_list_to_data_frame(self, list):
        """
        Convert lists with property to pandas dataframe
        """
        df = pd.DataFrame(list)
        return df

    def plot_overlap_density(
        self,
        property_name,
        df1,
        df2,
        outfilepath,
        n_breaks=50,
    ):

        min_plot_value = min(float(df1.min()), float(df2.min()))
        max_plot_value = max(float(df2.max()), float(df2.max()))

        breaks_seq = []
        for i in np.arange(
            min_plot_value, max_plot_value, (max_plot_value - min_plot_value) / n_breaks
        ):
            breaks_seq.append(i)

        plt.style.use("ggplot")
        plt.hist(df1, breaks_seq, alpha=0.5, label="reference", density=True)
        plt.hist(df2, breaks_seq, alpha=0.5, label="sample", density=True)
        plt.xlim(math.floor(min_plot_value), math.ceil(max_plot_value))
        plt.ylim(0, 1)
        plt.xlabel(
            "{}_{}".format(property_name, "/Å"),
            fontsize="12",
            horizontalalignment="center",
        )
        plt.ylabel(
            "density",
            fontsize="12",
            horizontalalignment="center",
        )
        plt.legend(loc="upper right")
        plt.savefig(outfilepath)

    def plot(self, prefix="property_vector", outfilepath=None):
        """
        Plots distribution of a given property vector and saves the file.

        Attributes
        -----------
        prefix: str
            Prefix that will be used in filename when saved in outfilepath
        outfilepath : str
            Path where output file with plot will be saved
        """
        if outfilepath is None:
            outfilepath = self.outfilepath

        reference = self.convert_list_to_data_frame(self.property.property_vector)
        sample = self.convert_list_to_data_frame(self.sampled_property.property_vector)

        self.plot_overlap_density(
            self.property_name, reference, sample, outfilepath, n_breaks=50
        )
