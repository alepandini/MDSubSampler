"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Namir Oues
    
    This file is part of the MDSubSampler software 
    (https://github.com/alepandini/MDSubSampler).
    Copyright (c) 2023 Namir Oues and Alessandro Pandini.

    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
    General Public License for more details.

    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
import random
from mdss.property import SampledProperty
from mdss.dissimilarity import *
from mdss.utilities import write_output_files
from mdss.utilities import plot_property


class ProteinSampler:
    """
    Class representing samplers used to sample a protein trajectory.

    Attributes
    ----------
    protein_property : ProteinProperty
        An instance of the ProteinProperty class representing the reference property.
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    output_folder : str
        Path to the output file where results will be saved, given as user input.
    file_prefix : str
        Prefix for output file naming, given as user input.
    dissimilarity_measure : Dissimilarity, optional
        An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.

    display_name : None
        Placeholder attribute. It is not used in the class implementation.
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
        """
        Initialize ProtienSampler object.

        Parameters
        ----------
        protein_property : ProteinProperty
            An instance of the ProteinProperty class representing the reference property.
        protein_data : ProteinData
            An instance of the ProteinData class representing the protein data.
        output_folder : str
            Path to the output file where results will be saved, given as user input.
        file_prefix : str
            Prefix for output file naming, given as user input.
        dissimilarity_measure : Dissimilarity, optional
            An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.
        """
        self.protein_property = protein_property
        self.protein_data = protein_data
        self.property_vector = protein_property.property_vector
        self.frame_indices = protein_data.frame_indices
        self.sampled_property_vector = None
        self.sampled_frame_indices = None
        self.samples_indices = []
        self.output_folder = output_folder
        self.file_prefix = file_prefix
        self.dissimilarity_measure = dissimilarity_measure

    def _create_data_list(self):
        """
        Create a list of data tuples from the property vector and frame indices.

        Returns
        -------
        list
            A list of data tuples where each tuple consists of an index and a property-index pair.
        """
        property_indices_tuples = list(zip(self.property_vector, self.frame_indices))
        data_list = [(i, v) for i, v in enumerate(property_indices_tuples)]
        return data_list

    def _create_sampled_property(self, sampled_data_vector, size):
        """
        Create a sampled property object based on the sampled data vector.

        Parameters
        ----------
        sampled_data_vector : list
            A list of sampled data tuples, where each tuple contains an index and a property-index pair.
        size : int
            The size of the sampled property.

        Returns
        -------
        SampledProperty
            A SampledProperty object representing the sampled property with the corresponding frame indices and sample indices.
        """
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
        Perform sampling for each size in the sample size list, and return the output.

        Parameters
        ----------
        perc_vector : list, optional
            A list of sample size percentages of the full trajectory size. Default is [25, 10, 5, 1, 0.5, 0.1].
        dissimilarity_threshold : int, optional
            The threshold for the dissimilarity measure to select the best sample size from the list. Default is None.
        step_recording : bool, optional
            Whether to return output files for all sample sizes from the list. Default is False.

        Returns
        -------
        SampledProperty or None
            The selected SampledProperty object that meets the dissimilarity threshold, or None if no sample meets the threshold.

        Notes
        -----
        If no sample meets the dissimilarity threshold, a warning message is printed, and None is returned.
        Otherwise, the selected SampledProperty object is returned.
        If percentage values provided are not less than 100, a warning message is printed.
        """

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
                        selected_sample_key = sampled_property.property_key
                if (
                    step_recording
                    or sampled_property.ref_dissimilarity <= dissimilarity_threshold
                ):
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
    Subclass of ProteinSampler class representing a random sampler for protein trajectories.

    Attributes
    ----------
    protein_property : ProteinProperty
        An instance of the ProteinProperty class representing the reference property.
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    output_folder : str
        Path to the output file where results will be saved, given as user input.
    file_prefix : str
        Prefix for output file naming, given as user input.
    seed_number : int, optional
        Seed number to initialize the random-number generator and ensure reproducibility. Default is 1999.
    dissimilarity_measure : Dissimilarity, optional
        An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.

    display_name : str
        The display name of the random sampler. Set to "Random Sampling".
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
        """
        Initialise RandomSampler object.

        Parameters
        ----------
        protein_property : ProteinProperty
            An instance of the ProteinProperty class representing the reference property.
        protein_data : ProteinData
            An instance of the ProteinData class representing the protein data.
        output_folder : str
            Path to the output file where results will be saved, given as user input.
        file_prefix : str
            Prefix for output file naming, given as user input.
        seed_number : int, optional
            Seed number to initialize the random-number generator and ensure reproducibility. Default is 1999.
        dissimilarity_measure : Dissimilarity, optional
            An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.]
        """
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
        Performs random sampling.

        Parameters
        ----------
        size : int
            Sample size specified by the user.

        Returns
        -------
        SampledProteinProperty
            An object of the SampledProteinProperty class containing the sampled data.
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
    Subclass of ProteinSampler class representing a stratified sampler for protein trajectories.

    Attributes
    ----------
    protein_property : ProteinProperty
        An instance of the ProteinProperty class representing the reference property.
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    output_folder : str
        Path to the output file where results will be saved, given as user input.
    file_prefix : str
        Prefix for output file naming, given as user input.
    strata_vector : list
        A 1D list containing the strata labels for each data point.
    dissimilarity_measure : Dissimilarity, optional
        An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.

    display_name : str
        The display name of the stratified sampler. Set to "Stratified Sampling".
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
        """
        Initialise StratifiedSampler object.

        Parameters
        ----------
        protein_property : ProteinProperty
            An instance of the ProteinProperty class representing the reference property.
        protein_data : ProteinData
            An instance of the ProteinData class representing the protein data.
        output_folder : str
            Path to the output file where results will be saved, given as user input.
        file_prefix : str
            Prefix for output file naming, given as user input.
        strata_vector : list
            A 1D list containing the strata labels for each data point.
        dissimilarity_measure : Dissimilarity, optional
            An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.
        """
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
        Performs stratified sampling.

        Parameters
        ----------
        size : int
            The sample size specified by the user.

        Returns
        -------
        SampledProperty
            A SampledProperty object representing the sampled protein property.

        Raises
        ------
        TypeError
            If the specified size is less than half the number of layers.

        Notes
        -----
        If strata size is smaller than required sample. Then a warnign message 'Sampling with replacement.'
        is printed.
        If strata vector is incosistent in size with property vector a warning message is printed.
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
    Subclass of ProteinSampler class representing a uniform sampler for protein trajectories.

    Attributes
    ----------
    protein_property : ProteinProperty
        An instance of the ProteinProperty class representing the reference property.
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    output_folder : str
        Path to the output file where results will be saved, given as user input.
    file_prefix : str
        Prefix for output file naming, given as user input.
    strata_number : int
        The number of intervals where sampling is done.
    dissimilarity_measure : Dissimilarity, optional
        An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.

    display_name : str
        The display name of the uniform sampler. Set to "Uniform Sampling".
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
        """
        Initialise UniformSampler object.

        Parameters
        ----------
        protein_property : ProteinProperty
            An instance of the ProteinProperty class representing the reference property.
        protein_data : ProteinData
            An instance of the ProteinData class representing the protein data.
        output_folder : str
            Path to the output file where results will be saved, given as user input.
        file_prefix : str
            Prefix for output file naming, given as user input.
        strata_number : int
            The number of intervals where sampling is done.
        dissimilarity_measure : Dissimilarity, optional
            An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.
        """
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
        Performs uniform sampling.

        Parameters
        ----------
        size : int
            The sample size specified by the user.

        Returns
        -------
        SampledProperty
            A SampledProperty object representing the sampled protein property.
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
    Subclass of ProteinSampler class representing a weighted sampler for protein trajectories.

    Attributes
    ----------
    protein_property : ProteinProperty
        An instance of the ProteinProperty class representing the reference property.
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    output_folder : str
        Path to the output file where results will be saved, given as user input.
    file_prefix : str
        Prefix for output file naming, given as user input.
    seed_number : int
        Seed number to initialize the random-number generator and ensure reproducibility. Default is 1999.
    weights_vector : list, optional
        Vector of weights for each element in the sample. If not provided, the weights will be estimated from the discretized property vector.
    dissimilarity_measure : Dissimilarity, optional
        An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.

    display_name : str
        The display name of the weighted sampler. Set to "Weighted Sampling".
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
        """
        Initialise WeightedSampler object.

        Parameters
        ----------
        protein_property : ProteinProperty
            An instance of the ProteinProperty class representing the reference property.
        protein_data : ProteinData
            An instance of the ProteinData class representing the protein data.
        output_folder : str
            Path to the output file where results will be saved, given as user input.
        file_prefix : str
            Prefix for output file naming, given as user input.
        seed_number : int
            Seed number to initialize the random-number generator and ensure reproducibility. Default is 1999.
        weights_vector : list, optional
            Vector of weights for each element in the sample. If not provided, the weights will be estimated from the discretized property vector.
        dissimilarity_measure : Dissimilarity, optional
            An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.
        """
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

            if not self.protein_property.discretized_property_vector:
                self.protein_property._property_statistics()
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
        Performs Weighted Sampling.

        Parameters
        ----------
        size : int
            The sample size specified by the user.

        Returns
        ----------
        SampledProperty
            A SampledProperty object representing the sampled protein property.

        Notes
        -----
        If weights vector has a different size from property vector, a warning message is printed.
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

            if len(data_list) == len(self.weights):
                sampled_data_vector = random.choices(
                    data_list, weights=self.weights, k=size
                )
                sampled_protein_property = self._create_sampled_property(
                    sampled_data_vector, size
                )
                return sampled_protein_property
            else:
                print("Error: The number of weights does not match the population.")


class BootstrappingSampler(ProteinSampler):
    """
    Subclass of ProteinSampler class representing a bootstrapping sampler for protein trajectories.

    Attributes
    ----------
    protein_property : ProteinProperty
        An instance of the ProteinProperty class representing the reference property.
    protein_data : ProteinData
        An instance of the ProteinData class representing the protein data.
    output_folder : str
        Path to the output file where results will be saved, given as user input.
    file_prefix : str
        Prefix for output file naming, given as user input.
    number_of_iterations : int
        Number of times the random sampling method is performed.
    seed_number : int
        Seed number to initialize the random-number generator and ensure reproducibility. Default is 1999.
    dissimilarity_measure : Dissimilarity, optional
        An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.

    display_name : str
        The display name of the bootstrapping sampler. Set to "Bootstrapping Sampling".
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
        """
        Initilise BootstrappingSampler object.

        Parameters
        ----------
        protein_property : ProteinProperty
            An instance of the ProteinProperty class representing the reference property.
        protein_data : ProteinData
            An instance of the ProteinData class representing the protein data.
        output_folder : str
            Path to the output file where results will be saved, given as user input.
        file_prefix : str
            Prefix for output file naming, given as user input.
        number_of_iterations : int
            Number of times the random sampling method is performed.
        seed_number : int
            Seed number to initialize the random-number generator and ensure reproducibility. Default is 1999.
        dissimilarity_measure : Dissimilarity, optional
            An instance of the Dissimilarity class representing the dissimilarity measure. Default is Bhattacharyya.
        """
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
        Performs Bootstrapping sampling.

        Parameters
        ----------
        size : int
            Sample size specified by the user.

        Returns
        -------
        SampledProteinProperty
            An object of the SampledProteinProperty class containing the sampled data.
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
