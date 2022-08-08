import numpy as np
import random


class ProteinSampler:
    """
    A class used to create a sample of a protein trajectory
    Attributes
    ----------
    protein_data: ProteinData class object
        The frame_list can be accessed through this object
    """

    display_name = None

    def __init__(self, property_vector):
        self.property_vector = list(property_vector)
        self.sampled_property_vector = None

    def sample(self, size):
        pass


class RandomSampler(ProteinSampler):
    """
     A Subclass of ProteinSampler class that uses Random Sampling
     Attributes
     ----------
    protein_data: ProteinData class object
         The frame_list can be accessed through this object
     seed: int
         Number that initialise a random-number generator
    """

    display_name = "Random Sampling"

    def __init__(self, property_vector, seed_number=1999):
        random.seed(seed_number)
        super().__init__(property_vector)

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
        temp_vector = [(i, v) for i, v in enumerate(self.property_vector)]
        temp_sampled_property_vector = random.sample(temp_vector, size)
        self.sampled_property_vector = [r[1] for r in temp_sampled_property_vector]
        self.sampled_indices = [r[0] for r in temp_sampled_property_vector]
        return self.sampled_property_vector


class UniformSampler(ProteinSampler):
    """
    A Subclass of ProteinSampler class that uses Uniform Sampling. A sample is generated
    from a frame list with uniform distribution.Samples are uniformly distributed over the
    half-open interval [low, high) (includes low, but excludes high)
    Attributes
    ----------
    protein_data: ProteinData class object
        The frame_list can be accessed through this object
    low: float
        Lower boundary of the output interval. The default value is 0.
    high: float
        Upper boundary of the output interval. The default value is 1.0.
    """

    display_name = "Uniform Sampling"

    def __init__(self, property_vector, low, high, dtype=int):
        self.low = low
        self.high = high
        self.dtype = dtype
        super().__init__(property_vector)

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
        self.sampled_property_vector = np.random.randint(
            self.low, self.high, size, self.dtype
        )
        super().sample(size)
        return self.sampled_property_vector


class StratifiedSampler(ProteinSampler):
    """
     A Subclass of ProteinSampler class that uses Stratified Sampling
    Attributes
         ----------
        protein_data: ProteinData class object
            The frame_list can be accessed through this object
         layers: list
            2D list that consists of multiple layers where eachlayer is a set of
            labels for the frames according to the strata
    """

    display_name = "Stratified Sampling"

    def __init__(self, property_vector, layers):
        self.layers = layers
        super().__init__(property_vector)

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
    protein_data: object
        ProteinData class object that has access to all methods and attributes
        of ProteinData class. The frame_list can be accessed through it.
    number_of_iterations: int
        This is the number of times the random sampling method is performed
    """

    display_name = "Bootstrapping Sampling"

    def __init__(self, property_vector, number_of_iterations):
        self.number_of_iterations = number_of_iterations
        super().__init__(property_vector)

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
            current_sample = np.random.choice(self.property_vector, size, replace=True)
            samples.append(current_sample)
        return samples
