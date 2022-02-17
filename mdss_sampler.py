import numpy as np
import random


class ProteinSampler:
    """
    A class used to create a sample of a protein trajectory

    Attributes
    ----------
    frame_list: list
        List that contains all the frames from a given protein trajectory
    """

    def __init__(self, frame_list):
        self.frame_list = frame_list
        self.sampled_frame_list = None

    def sample(self, size, number_of_iterations=None):
        self.sampled_frame_list.sort()


class RandomSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Random Sampling

    Attributes
    ----------
    frame_list: list
        List that contains all the frames from a given protein trajectory
    seed: int
        Number that initialise a random-number generator
    """

    def __init__(self, frame_list, seed_number=1999):
        random.seed(seed_number)
        super().__init__(frame_list)

    def sample(self, size, number_of_iterations=None):
        """
        Method that generates a random sample of a list

        Attributes
        ----------
        size: int
            The sample size

        Returns
        ----------
        return a single random sample of the frame list with the desired size
        """
        self.sampled_frame_list = random.sample(self.frame_list, size)
        super().sample(size)
        return self.sampled_frame_list


class StratifiedSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Stratified Sampling
    """

    def __init__(self, frame_list, layers):
        self.layers = layers
        super().__init__(frame_list)

    def strata_sample_size(self, size, population, layer_size):
        """
        Method that calculates the sample size of the strata (ie homogeneous groups)

        Attributes
        ----------
        size: int
            Whole sample size
        population: int
            Whole population size (sum of each layer - ie strata - size)
        layer_size: int
            Size for current layer

        Returns
        ----------
        return the rounded sample size of the strata

        """
        cur_size = (size / population) * layer_size
        return round(cur_size)

    def sample(self, size, number_of_iterations=None):
        """
        Method that does the stratified sampling

        Attributes
        ----------
        layers: vector
            This is a 2D vector that consists of multiple layers. Each layer is a set of labels
            for the frames according to the strata
        size: int
            Whole sample size

        Returns
        ----------
        return a single stratified sample

        """
        population = sum(len(layer) for layer in self.layers)
        samples = []

        for layer in self.layers:
            layer_size = len(layer)
            current_layer_sample_size = self.strata_sample_size(
                size, population, layer_size
            )
            current_sample = random.sample(layer, current_layer_sample_size)
            samples.extend(current_sample)
        return samples


class BootstrappingSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Bootstrapping Sampling

    Attributes
    ----------
    frame_list: list
        List that contains all the frames from a given protein trajectory
    """

    def sample(self, size, number_of_iterations):
        """
        Method that does the bootstrapping sampling

        Attributes
        ----------
        size: int
            This is the desired size of the sample each time we iterate
        number_of_iterations: int
            This is the number of times the random sampling method is performed

        Returns
        ----------
        return a list of bootstrapped samples

        """
        samples = []
        for i in range(number_of_iterations):
            current_sample = np.random.choice(self.frame_list, size, replace=True)
            samples.append(current_sample)
        return samples


# samples = []
#         for i in range(number_of_iterations):
#             current_sample = np.random.choice(self.frame_list, size, replace=True)
#             current_sample_mean = self.find_nearest(current_sample, np.mean(current_sample))
#             samples.append(current_sample_mean)

#         return samples

# def find_nearest(self, array, value):
#     """
#     Method that finds the closest number to the mean value

#     Attributes
#     ----------
#     array: list
#         This is a list with the data
#     value: float
#         This is the number of times the random sampling method is performed

#     Returns
#     ----------
#     return an integer with the frame number that is closest to the mean number

#     """
#     idx = (np.abs(array - value)).argmin()
#     return array[idx]
