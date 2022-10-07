import numpy as np
import random
import mdss_property
from mdss_property import SampledProperty


class ProteinSampler:
    """
    A class used to create a sample of a protein trajectory
    Attributes
    ----------
    protein_data: ProteinData class object
        The frame_list can be accessed through this object
    """

    display_name = None

    def __init__(self, protein_property):
        self.protein_property = protein_property
        self.property_vector = protein_property.property_vector
        self.frame_indices = protein_property.frame_indices
        self.sampled_property_vector = None
        self.sampled_frame_indices = None

    def _create_data_list(self):
        property_indices_tuples = list(zip(self.property_vector, self.frame_indices))
        data_list = [(i, v) for i, v in enumerate(property_indices_tuples)]
        return data_list

    def _create_sampled_property(self, sampled_data_vector):
        self.sampled_property_vector = [r[1][0] for r in sampled_data_vector]
        self.sampled_frame_indices = [r[1][1] for r in sampled_data_vector]
        sampled_protein_property = SampledProperty(
            self.protein_property,
            self.sampled_property_vector,
            self.sampled_frame_indices,
        )
        return sampled_protein_property

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

    def __init__(self, protein_property, seed_number=1999):
        random.seed(seed_number)
        super().__init__(protein_property)

    def sample(self, size):
        """
        Method that generates a random sample of a list
        Attributes
        ----------
        size: int
            The sample size ?
        Returns
        ----------
        return a single random sample of the frame list with the desired size
        """
        data_list = self._create_data_list()
        sampled_data_vector = random.sample(data_list, size)
        sampled_protein_property = self._create_sampled_property(sampled_data_vector)
        return sampled_protein_property


class StratifiedSampler(ProteinSampler):
    """
     A Subclass of ProteinSampler class that uses Stratified Sampling
    Attributes
         ----------
        protein_data: ProteinData class object
            The frame_list can be accessed through this object
         layers: list
            2D list that consists of strata lables one for each point
    """

    display_name = "Stratified Sampling"

    def __init__(self, protein_property, strata_vector):
        self.strata_vector = strata_vector
        strata_labels = sorted(set(strata_vector))
        self.layers = {}
        for label in strata_labels:
            strata_indices = [
                idx for idx, value in enumerate(strata_vector) if value == label
            ]
            self.layers[label] = strata_indices
        self.n_layers = len(self.layers.keys())
        super().__init__(protein_property)

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
        population_size = sum(len(layer) for layer in self.layers.values())
        if population_size == len(self.protein_property.property_vector):
            sampled_data_vector = []
            strata_sample_size = round(size / self.n_layers)
            if strata_sample_size < 1:
                print("Warning: size should be at least half the number of layers.")
                return None

            for layer in self.layers.values():
                if len(layer) < strata_sample_size:
                    print(
                        "Warning: strata size smaller than required sample. Sampling with replacement."
                    )
                    layer_sampled_indices = random.choices(layer, k=strata_sample_size)
                else:
                    layer_sampled_indices = random.sample(layer, strata_sample_size)
                sampled_property_vector = [
                    self.property_vector[x] for x in layer_sampled_indices
                ]
                sampled_frame_indices = [
                    self.frame_indices[x] for x in layer_sampled_indices
                ]
                property_indices_tuples = list(
                    zip(sampled_property_vector, sampled_frame_indices)
                )
                sampled_data_vector.extend(
                    list(zip(layer_sampled_indices, property_indices_tuples))
                )

            sampled_protein_property = self._create_sampled_property(
                sampled_data_vector
            )
            return sampled_protein_property
        else:
            print("Warning: strata vector is incosistent in size with property vector.")
            return None


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

    def __init__(self, protein_property, strata_number):
        super().__init__(protein_property)
        self.bin_size = (
            self.protein_property.max_value - self.protein_property.min_value
        ) / strata_number
        self.bins_vector = np.arange(
            self.protein_property.min_value,
            self.protein_property.max_value,
            self.bin_size,
        )
        self.strata_vector = np.digitize(self.property_vector, self.bins_vector)

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
        strat_sampler = StratifiedSampler(self.protein_property, self.strata_vector)
        sampled_protein_property = strat_sampler.sample(size)
        return sampled_protein_property


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

    def __init__(self, protein_property, number_of_iterations, seed_number=1999):
        random.seed(seed_number)
        self.number_of_iterations = number_of_iterations
        super().__init__(protein_property)

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
        data_list = self._create_data_list()
        sampled_data_vector = []
        for i in range(self.number_of_iterations):
            current_sample = random.sample(data_list, size)
            sampled_data_vector.extend(current_sample)
        sampled_protein_property = self._create_sampled_property(sampled_data_vector)
        return sampled_protein_property
