import os


class PropertyVectorPlot:
    """
    Represents the vector plot after property calculation

    Attributes
    -----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    outfilepath : str
        path where output file with plot will be saved
    """

    def __init__(self, protein_property, outfilepath):
        self.protein_property = protein_property
        self.property_vector = protein_property.property_vector
        self.frame_indices = protein_property.frame_indices
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
        pass
