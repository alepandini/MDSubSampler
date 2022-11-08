import os


class PropertyVectorPlot:

    def __init__(self, protein_property, outfilepath):
        self.protein_property = protein_property
        self.property_vector = protein_property.property_vector
        self.frame_indices = protein_property.frame_indices
        self.outfilepath = outfilepath

    def plot(self, prefix="property_vector", outfilepath=None):
        if outfilepath is None:
            outfilepath = self.outfilepath
        pass
