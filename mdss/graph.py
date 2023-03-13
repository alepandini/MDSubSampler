import os
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import MaxAbsScaler
from matplotlib import style


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
        self.sampled_property_vector = sampled_property.property_vector
        self.frame_indices = property.frame_indices
        self.outfilepath = outfilepath

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
        # Convert lists with property to pandas dataframe
        sample = self.sampled_property.property_vector
        reference = self.property.property_vector
        sample = pd.DataFrame(sample)
        reference = pd.DataFrame(reference)

        # Apply scaling to dataframes
        max = MaxAbsScaler()
        sample = max.fit_transform(sample)
        reference = max.fit_transform(reference)

        # Plot the density distribution of the data
        plt.style.use("ggplot")
        plt.hist(reference, bins=20, alpha=0.5, density=True, label="reference")
        plt.hist(sample, bins=20, alpha=0.5, density=True, label="sample")
        plt.legend(loc="upper right")
        plt.savefig(outfilepath)
