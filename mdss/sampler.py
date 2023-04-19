import numpy as np
import random
from mdss.property import SampledProperty
from mdss.dissimilarity import *
from mdss.utilities import write_output_files
from mdss.utilities import plot_property


class ProteinSampler:
    """
    Represents Sampler Class

    Attributes
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    output_folder: str
        Path to output file given as user input
    file_prefix: str,
        Prefix given as user input
    dissimilarity_measure: Dissimilarity class object
        Dissimilarity measure between full and sample trajectory
    """

    display_name = None

    def __init__(
        self,
        protein_property,
        protein_data,
        output_folder,
        file_prefix,
        dissimilarity_measure=Bhattacharyya,
    ):
        self.protein_property = protein_property
        self.protein_data = protein_data
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

    def _create_sampled_property(self, sampled_data_vector, size):
        self.sampled_property_vector = [r[1][0] for r in sampled_data_vector]
        self.sampled_frame_indices = [r[1][1] for r in sampled_data_vector]
        sampled_protein_property = SampledProperty(
            self.protein_property,
            self.sampled_property_vector,
            self.sampled_frame_indices,
            self.samples_indices,
            size,
            self.dissimilarity_measure,
        )
        return sampled_protein_property

    def sample(self, size):
        pass

    def scan_sample_size(
        self, perc_vector=None, dissimilarity_threshold=None, step_recording=False
    ):
        """
         Scans sample size list (user input), performs sampling for each size and returns output

         Attributes
         ----------
         perc_vector: list
             Contains list of sample size percentages of full trajectory size
         dissimilarity_threshold: int
             Threshold for dissimilarity measure to get the best sample size from the list
        step_recording: Boolean
             Returns output files for all sample sizes from the list
        """
        selected_size = None
        selected_sample_key = None
        if perc_vector is None:
            perc_vector = [25, 10, 5, 1, 0.5, 0.1]
        if sum([x > 100 for x in perc_vector]) == 0:
            n_frames = len(self.property_vector)
            perc_vector.sort(reverse=True)
            for p in perc_vector:
                log.info(
                    "{:15s} Sample percentage size: {:4.5s}".format(
                        "INPUT", str(round(p))
                    )
                )
                sampled_property = self.sample(round(p * n_frames / 100))
                if sampled_property is not None:
                    if dissimilarity_threshold is None:
                        dissimilarity_threshold = (
                            sampled_property.dissimilarity_threshold
                        )
                    if sampled_property.ref_dissimilarity <= dissimilarity_threshold:
                        selected_size = p
                        selected_sample_key = sampled_property.property_key
                if step_recording:
                    write_output_files(
                        output_folder=self.output_folder,
                        file_prefix=self.file_prefix,
                        p_prop=self.protein_property,
                        s_prop=sampled_property,
                        p_data=self.protein_data,
                        p=p,
                    )
                    plot_property(
                        output_folder=self.output_folder,
                        file_prefix=self.file_prefix,
                        p_prop=self.protein_property,
                        s_prop=sampled_property,
                        p=p,
                    )
                else:
                    write_output_files(
                        output_folder=self.output_folder,
                        file_prefix=self.file_prefix,
                        p_prop=self.protein_property,
                        s_prop=sampled_property,
                        p_data=self.protein_data,
                    )
                    plot_property(
                        output_folder=self.output_folder,
                        file_prefix=self.file_prefix,
                        p_prop=self.protein_property,
                        s_prop=sampled_property,
                    )

            log.info(
                "{:15s} Output files for all sample sizes were generated successfully".format(
                    "OUTPUT"
                )
            )
            if selected_sample_key is None:
                print("Warning: no sample meeting dissimilarity threshold")
                log.warning(
                    "{:12s} No sample meeting dissimilarity threshold".format("STEPS")
                )
                return None
            else:
                return self.protein_property.protein_data.property_dict[
                    selected_sample_key
                ]
        else:
            print("Percentage values should be smaller than 100")
            log.warning(
                "{:12s} Percentage values should be smaller than 100".format("INPUT")
            )


class RandomSampler(ProteinSampler):
    """
    Represents Random Sampler

    Attributes
     ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    output_folder: str
        Path to output file given as user input
    file_prefix: str,
        Prefix given as user input
    seed_number: int,
       Initialise a random-number generator and ensures reproducibility
    dissimilarity_measure: Dissimilarity class object
        Dissimilarity measure between full and sample trajectory
    """

    display_name = "Random Sampling"

    def __init__(
        self,
        protein_property,
        protein_data,
        output_folder,
        file_prefix,
        seed_number=1999,
        dissimilarity_measure=Bhattacharyya,
    ):
        random.seed(seed_number)
        super().__init__(
            protein_property=protein_property,
            protein_data=protein_data,
            output_folder=output_folder,
            file_prefix=file_prefix,
            dissimilarity_measure=dissimilarity_measure,
        )

    def sample(self, size):
        """
        Performs Random Sampling
        Attributes
        ----------
        size: int
            Sample size - user input

        Returns
        ----------
        A sampled_protein_property object
        """
        self.samples_indices = list(np.repeat(0, size))
        data_list = self._create_data_list()
        sampled_data_vector = random.sample(data_list, size)
        sampled_protein_property = self._create_sampled_property(
            sampled_data_vector, size
        )
        return sampled_protein_property


class StratifiedSampler(ProteinSampler):
    """
    Represents Stratified sampler class
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    output_folder: str
        Path to output file given as user input
    file_prefix: str,
        Prefix given as user input
    strata_vector: list
        2D list that consists of strata labels one for each point
    dissimilarity_measure: Dissimilarity class object
        Dissimilarity measure between full and sample trajectory
    """

    display_name = "Stratified Sampling"

    def __init__(
        self,
        protein_property,
        protein_data,
        output_folder,
        file_prefix,
        strata_vector,
        dissimilarity_measure=Bhattacharyya,
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
            protein_property=protein_property,
            protein_data=protein_data,
            output_folder=output_folder,
            file_prefix=file_prefix,
            dissimilarity_measure=dissimilarity_measure,
        )

    def sample(self, size):
        """
        Performs Stratified Sampling
        Attributes
        ----------
        size: int
            Sample size user input

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
                log.warn(
                    "{:12s} Size should be at least half the number of layers".format(
                        "INPUT"
                    )
                )
                raise TypeError("Size should be at least half the number of layers")

            self.samples_indices = []
            sample_index = 0
            for layer in self.layers.values():
                if len(layer) < strata_sample_size:
                    print(
                        "Warning: strata size smaller than required sample. Sampling with replacement."
                    )
                    log.warning(
                        "{:12s} Strata size smaller than required sample. Sampling with replacement.".format(
                            "INPUT"
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
                sampled_data_vector, size
            )
            return sampled_protein_property
        else:
            print("Warning: strata vector is incosistent in size with property vector")
            log.warning(
                "{:12s} Strata vector is incosistent in size with property vector".format(
                    "STEPS"
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
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    output_folder: str
        Path to output file given as user input
    file_prefix: str,
        Prefix given as user input
    strata_number: int
        The number of intervals where sampling is done.
    dissimilarity_measure: Dissimilarity class object
        Dissimilarity measure between full and sample trajectory
    """

    display_name = "Uniform Sampling"

    def __init__(
        self,
        protein_property,
        protein_data,
        output_folder,
        file_prefix,
        strata_number,
        dissimilarity_measure=Bhattacharyya,
    ):
        super().__init__(
            protein_property=protein_property,
            protein_data=protein_data,
            output_folder=output_folder,
            file_prefix=file_prefix,
            dissimilarity_measure=dissimilarity_measure,
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

    def sample(self, size):
        """
        Performs Uniform Sampling

        Returns
        ----------
        A sampled_protein_property object
        """
        strat_sampler = StratifiedSampler(
            self.protein_property,
            self.protein_data,
            self.output_folder,
            self.file_prefix,
            self.strata_vector,
            self.dissimilarity_measure,
        )
        sampled_protein_property = strat_sampler.sample(size)
        return sampled_protein_property


class WeightedSampler(ProteinSampler):
    """
    Represents Wighted Sampler
    ----------
    protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    output_folder: str
        Path to output file given as user input
    file_prefix: str,
        Prefix given as user input
    seed_number: int,
       Initialise a random-number generator and ensures reproducibility
    weights_vector:
        Vector of weights for each element in the sample
    dissimilarity_measure: Dissimilarity class object
        Dissimilarity measure between full and sample trajectory
    """

    display_name = "Weighted Sampling"

    def __init__(
        self,
        protein_property,
        protein_data,
        output_folder,
        file_prefix,
        seed_number=1999,
        weights_vector=None,
        dissimilarity_measure=Bhattacharyya,
    ):
        random.seed(seed_number)
        super().__init__(
            protein_property=protein_property,
            protein_data=protein_data,
            output_folder=output_folder,
            file_prefix=file_prefix,
            dissimilarity_measure=dissimilarity_measure,
        )
        if weights_vector is None:
            print(
                "Weights not provided. They will be estimated from discretized property vector."
            )
            log.warning(
                "{:12s} Weights not provided. They will be estimated from discretized property vector.".format(
                    "STEPS"
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

    def sample(self, size):
        """
        Performs Weighted Sampling

        Attributes
        ----------
        size: int
            Sample size - user input

        Returns
        ----------
        A sampled_protein_property object
        """
        if len(self.weights) != len(self.property_vector):
            print("Warning: weights vector of different size from property vector")
            log.warning(
                "{:12s} Weights vector of different size from property vector".format(
                    "STEPS"
                )
            )
        else:
            self.samples_indices = list(np.repeat(0, size))
            data_list = self._create_data_list()
            sampled_data_vector = random.choices(
                data_list, weights=self.weights, k=size
            )
            sampled_protein_property = self._create_sampled_property(
                sampled_data_vector, size
            )
            return sampled_protein_property


class BootstrappingSampler(ProteinSampler):
    """
    Represents Bootstrapping Sampler

    Attributes
    ----------
     protein_property: ProteinProperty class object
        The object has access to all methods and attributes of ProteinProperty class
    protein_data: ProteinData class object
        The object has access to all methods and attributes of ProteinData class
    output_folder: str
        Path to output file given as user input
    file_prefix: str,
        Prefix given as user input
    number_of_iterations: int
        Number of times the random sampling method is performed
    seed_number: int,
       Initialise a random-number generator and ensures reproducibility
    dissimilarity_measure: Dissimilarity class object
        Dissimilarity measure between full and sample trajectory
    """

    display_name = "Bootstrapping Sampling"

    def __init__(
        self,
        protein_property,
        protein_data,
        output_folder,
        file_prefix,
        number_of_iterations,
        seed_number=1999,
        dissimilarity_measure=Bhattacharyya,
    ):
        random.seed(seed_number)
        self.number_of_iterations = number_of_iterations
        super().__init__(
            protein_property=protein_property,
            protein_data=protein_data,
            output_folder=output_folder,
            file_prefix=file_prefix,
            dissimilarity_measure=dissimilarity_measure,
        )

    def sample(self, size):
        """
        Performs Bootstrapping sampling

        Attributes
        ----------
        size: int
            Sample size - user input

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
        sampled_protein_property = self._create_sampled_property(
            sampled_data_vector, size
        )
        return sampled_protein_property
