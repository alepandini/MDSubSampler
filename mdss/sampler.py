import numpy as np
import random
from mdss.property import SampledProperty
from mdss.dissimilarity import *

import sys
import os


def convert_size(size, n_frames):
    """
    Converts sample size of traj into int
    User input could be given as percentage or int
    """
    if isinstance(size, int):
        return size

    if size.endswith("%"):
        prc = float(size.rstrip("%"))
        if prc > 100:
            print("size percentage {} is not less than 100%".format(prc))
            log.error(
                "{:15s} Size percentage {} is not less than 100%".format("ERROR", prc)
            )
            sys.exit(1)

        size = prc * n_frames / 100
        return round(size)

    log.info("{:15s} Sample size was converted to int".format("STEPS"))
    return int(size)


class ProteinSampler:
    """
    Represents Sampler Class

    Attributes
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    dissimilarity_measure: Dissimilarity class object
        Dissimilarity measure between full and sample trajectory
    """

    display_name = None

    def __init__(
        self,
        protein_property,
        output_folder,
        file_prefix,
        dissimilarity_measure=Bhattacharya,
    ):
        self.protein_property = protein_property
        self.property_vector = protein_property.property_vector
        self.frame_indices = protein_property.frame_indices
        self.sampled_property_vector = None
        self.sampled_frame_indices = None
        self.samples_indices = []
        self.output_folder = output_folder
        self.file_prefix = file_prefix
        self.dissimilarity_measure = dissimilarity_measure

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
            self.samples_indices,
            self.dissimilarity_measure,
        )
        log.info("{:15s} Sampled property was created".format("STEPS"))
        return sampled_protein_property

    def sample(self, size):
        return self._sample(
            convert_size(size, self.protein_property.protein_data.n_frames)
        )

    def _sample(self, size):
        pass

    def scan_sample_size(
        self, perc_vector=None, dissimilarity_threshold=None, step_recording=True
    ):
        selected_size = None
        selected_sample_key = None
        if perc_vector is None:
            perc_vector = [25, 10, 5, 1, 0.5, 0.1]
        if sum([x > 100 for x in perc_vector]) == 0:
            n_frames = len(self.property_vector)
            perc_vector.sort(reverse=True)
            if step_recording == True:
                for p in perc_vector:
                    sampled_property = self.sample(round(p * n_frames / 100))
                    if sampled_property is not None:
                        if dissimilarity_threshold is None:
                            dissimilarity_threshold = (
                                sampled_property.dissimilarity_threshold
                            )
                        if (
                            sampled_property.ref_dissimilarity
                            <= dissimilarity_threshold
                        ):
                            selected_size = p
                            selected_sample_key = sampled_property.property_key
                    filename = "{}_{}_{}_{}.dat".format(
                        self.file_prefix,
                        p,
                        self.protein_property.display_name,
                        self.dissimilarity_measure.display_name,
                    )
                    filepath = os.path.join(self.output_folder, filename)
                    sampled_property.write_property_vector(filepath)
            if selected_sample_key is None:
                print("Warning: no sample meeting dissimilarity threshold")
                log.warning(
                    "{:15s} No sample meeting dissimilarity threshold".format("WARNING")
                )
                return None
            else:
                return self.protein_property.protein_data.property_dict[
                    selected_sample_key
                ]
        else:
            print("Percentage values should be smaller than 100")
            log.warning(
                "{:15s} Percentage values should be smaller than 100".format("WARN")
            )


class RandomSampler(ProteinSampler):
    """
    Represents Random Sampler

    Attributes
     ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    seed: int
        Number that initialise a random-number generator
    dissimilarity_measure: Dissimilarity class object
        Default measure is Bhattacharya
    """

    display_name = "Random Sampling"

    def __init__(
        self,
        protein_property,
        output_folder,
        file_prefix,
        seed_number=1999,
        dissimilarity_measure=Bhattacharya,
    ):
        random.seed(seed_number)
        super().__init__(
            protein_property,
            output_folder,
            file_prefix,
            dissimilarity_measure,
        )

    def _sample(self, size):
        """
        Performs Random Sampling
        Attributes
        ----------
        size: int
            The sample size - user's choice
        Returns
        ----------
        A sampled_protein_property object
        """
        self.samples_indices = list(np.repeat(0, size))
        data_list = self._create_data_list()
        sampled_data_vector = random.sample(data_list, size)
        sampled_protein_property = self._create_sampled_property(sampled_data_vector)
        log.info("{:15s} Random sample was created".format("STEPS"))
        return sampled_protein_property


class StratifiedSampler(ProteinSampler):
    """
    Represents Stratified sampler class
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    strata_vector: list
        2D list that consists of strata lables one for each point
    dissimilarity_measure: Dissimilarity class object
        Default measure is Bhattacharya
    """

    display_name = "Stratified Sampling"

    def __init__(
        self,
        protein_property,
        output_folder,
        file_prefix,
        strata_vector,
        dissimilarity_measure=Bhattacharya,
    ):
        self.strata_vector = strata_vector
        strata_labels = sorted(set(strata_vector))
        self.layers = {}
        for label in strata_labels:
            strata_indices = [
                idx for idx, value in enumerate(strata_vector) if value == label
            ]
            self.layers[label] = strata_indices
        self.n_layers = len(self.layers.keys())
        super().__init__(
            protein_property,
            output_folder,
            file_prefix,
            dissimilarity_measure,
        )

    def _sample(self, size):
        """
        Performs Stratified Sampling
        Attributes
        ----------
        size: int
            The sample size - user's choice
        Returns
        ----------
        A sampled_protein_property object
        """
        population_size = sum(len(layer) for layer in self.layers.values())
        if population_size == len(self.protein_property.property_vector):
            sampled_data_vector = []
            strata_sample_size = round(size / self.n_layers)
            if strata_sample_size < 1:
                print("Warning: Size should be at least half the number of layers")
                log.warning(
                    "{:15s} Size should be at least half the number of layers".format(
                        "WARNING"
                    )
                )
                return None

            self.samples_indices = []
            sample_index = 0
            for layer in self.layers.values():
                if len(layer) < strata_sample_size:
                    print(
                        "Warning: strata size smaller than required sample. Sampling with replacement."
                    )
                    log.warning(
                        "{:15s} Strata size smaller than required sample. Sampling with replacement.".format(
                            "WARNING"
                        )
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
                self.samples_indices.extend(
                    list(np.repeat(sample_index, strata_sample_size))
                )
                sample_index += 1

            sampled_protein_property = self._create_sampled_property(
                sampled_data_vector
            )
            log.info("{:15s} Stratified sample was created".format("STEPS"))
            return sampled_protein_property
        else:
            print("Warning: strata vector is incosistent in size with property vector")
            log.warning(
                "{:15s} Strata vector is incosistent in size with property vector".format(
                    "WARNING"
                )
            )
            return None


class UniformSampler(ProteinSampler):
    """
    Represents Uniform Sampler
    Attributes
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    strata_number: int
        The number of intervals where sampling is done.
    dissimilarity_measure: Dissimilarity class object
        Default measure is Bhattacharya
    """

    display_name = "Uniform Sampling"

    def __init__(
        self,
        protein_property,
        output_folder,
        file_prefix,
        strata_number,
        dissimilarity_measure=Bhattacharya,
    ):
        super().__init__(
            protein_property,
            output_folder,
            file_prefix,
            dissimilarity_measure,
        )
        self.bin_size = (
            self.protein_property.max_value - self.protein_property.min_value
        ) / strata_number
        self.bins_vector = np.arange(
            self.protein_property.min_value,
            self.protein_property.max_value,
            self.bin_size,
        )
        self.strata_vector = np.digitize(self.property_vector, self.bins_vector)

    def _sample(self, size):
        """
        Performs Uniform Sampling
        Returns
        ----------
        A sampled_protein_property object
        """
        strat_sampler = StratifiedSampler(self.protein_property, self.strata_vector)
        sampled_protein_property = strat_sampler.sample(size)
        log.info("{:15s} Uniform sample was created".format("STEPS"))
        return sampled_protein_property


class WeightedSampler(ProteinSampler):
    """
    Represents Wighted Sampler
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    seed: int
        Number that initialise a random-number generator
    weights_vector:
        Vector of weights for each element in the sample
    dissimilarity_measure: Dissimilarity class object
        Default measure is Bhattacharya
    """

    display_name = "Weighted Sampling"

    def __init__(
        self,
        protein_property,
        output_folder,
        file_prefix,
        seed_number=1999,
        weights_vector=None,
        dissimilarity_measure=Bhattacharya,
    ):
        random.seed(seed_number)
        super().__init__(
            protein_property,
            output_folder,
            file_prefix,
            dissimilarity_measure,
        )
        if weights_vector is None:
            print(
                "Weights not provided. They will be estimated from discretized property vector."
            )
            log.warning(
                "{:15s} Weights not provided. They will be estimated from discretized property vector.".format(
                    "WARNING"
                )
            )
            if self.protein_property.discretized_property_vector is None:
                self.protein_property.discretize_vector()
            self.weights = []
            for value in self.protein_property.discretized_property_vector:
                self.weights.append(
                    self.protein_property.discretized_property_vector.count(value)
                )
        else:
            self.weights = weights_vector

    def _sample(self, size):
        """
        Performs Weighted Sampling
        Attributes
        ----------
        size: int
            The sample size - user's choice
        Returns
        ----------
        A sampled_protein_property object
        """
        if len(self.weights) != len(self.property_vector):
            print("Warning: weights vector of different size from property vector")
            log.warning(
                "{:15s} Weights vector of different size from property vector".format(
                    "WARNING"
                )
            )
        else:
            self.samples_indices = list(np.repeat(0, size))
            data_list = self._create_data_list()
            sampled_data_vector = random.choices(
                data_list, weights=self.weights, k=size
            )
            sampled_protein_property = self._create_sampled_property(
                sampled_data_vector
            )
            log.info("{:15s} Weighted sample was created".format("STEPS"))
            return sampled_protein_property


class BootstrappingSampler(ProteinSampler):
    """
    Represents Bootstrapping Sampler
    Attributes
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    number_of_iterations: int
        This is the number of times the random sampling method is performed
    seed: int
        Number that initialise a random-number generator
    dissimilarity_measure: Dissimilarity class object
        Default measure is Bhattacharya
    """

    display_name = "Bootstrapping Sampling"

    def __init__(
        self,
        protein_property,
        output_folder,
        file_prefix,
        number_of_iterations,
        seed_number=1999,
        dissimilarity_measure=Bhattacharya,
    ):
        random.seed(seed_number)
        self.number_of_iterations = number_of_iterations
        super().__init__(
            protein_property,
            output_folder,
            file_prefix,
            dissimilarity_measure,
        )

    def _sample(self, size):
        """
        Performs Bootstrapping sampling
        Attributes
        ----------
        size: int
            The sample size - user's choice
        Returns
        ----------
        A sampled_protein_property object
        """
        data_list = self._create_data_list()
        sampled_data_vector = []
        self.samples_indices = []
        for i in range(self.number_of_iterations):
            current_sample = random.sample(data_list, size)
            sampled_data_vector.extend(current_sample)
            self.samples_indices.extend(list(np.repeat(i, size)))
        sampled_protein_property = self._create_sampled_property(sampled_data_vector)
        log.info("{:15s} Bootstrapped sample was created".format("STEPS"))
        return sampled_protein_property
