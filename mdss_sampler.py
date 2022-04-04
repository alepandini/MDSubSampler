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

    def sample(self, size):
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

    def sample(self, size):
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


class UniformSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Uniform Sampling. A sample is generated
    from a frame list with uniform distribution.Samples are uniformly distributed over the
    half-open interval [low, high) (includes low, but excludes high)

    Attributes
    ----------
    frame_list: list
        List that contains all the frames from a given protein trajectory
    low: float
        Lower boundary of the output interval. The default value is 0.

    high: float
        Upper boundary of the output interval. The default value is 1.0.

    """

    def __init__(self, frame_list, low, high, dtype=int):
        self.low = low
        self.high = high
        self.dtype = dtype
        super().__init__(frame_list)

    def sample(self, size):
        """
        Method that generates a uniform sample of a list

        Attributes
        ----------
        size: int
            The sample size

        Returns
        ----------
        Return random integers from the “discrete uniform” distribution of the specified dtype
        in the “half-open” interval [low, high).
        """
        self.sampled_frame_list = np.random.randint(
            self.low, self.high, size, self.dtype
        )
        super().sample(size)
        return self.sampled_frame_list


class StratifiedSampler(ProteinSampler):
    """
     A Subclass of ProteinSampler class that uses Stratified Sampling

    Attributes
         ----------
         frame_list: int
             List that contains all the frames from a given protein trajectory
         layers: list
            2D list that consists of multiple layers where eachlayer is a set of
            labels for the frames according to the strata
    """

    def __init__(self, frame_list, layers):
        self.layers = layers
        super().__init__(frame_list)

    def sample(self, size):
        """
        Method that performs the stratified sampling

        Attributes
        ----------
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
            # the sample size of the strata (ie homogeneous groups)
            cur_size = round((size / population) * layer_size)
            current_layer_sample_size = cur_size
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
    number_of_iterations: int
        This is the number of times the random sampling method is performed
    """

    def __init__(self, frame_list, number_of_iterations):
        self.number_of_iterations = number_of_iterations
        super().__init__(frame_list)

    def sample(self, size):
        """
        Method that does the bootstrapping sampling

        Attributes
        ----------
        size: int
            This is the desired size of the sample each time we iterate

        Returns
        ----------
        return a list of bootstrapped samples

        """
        samples = []
        for i in range(self.number_of_iterations):
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
