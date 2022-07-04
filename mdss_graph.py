import os


def plot_distribution(property, data_filepath, prefix, output_folder):
    filename = "{}_plot.png".format(prefix)
    filepath = os.path.join(output_folder, filename)
    property.plot_property(data_filepath, filepath)
