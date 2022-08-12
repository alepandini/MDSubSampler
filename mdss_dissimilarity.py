import dictances


class Dissimilarity:
    """
    A class used to calculate the dissimilarity distance in terms of property between
    a full protein trajectory and a sample of it to identify the difference

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = None

    def __init__(self, property_1, property_2, clean=False):
        self.property_1 = property_1
        self.property_2 = property_2
        self.dissimilarity = self.calculate_dissimilarity()

    def calculate_dissimilarity(self):
        """
        Method that calculates the difference between the average values of the
        two calculated property vectors.
        """
        return self.property_1.avg_value - self.property_2.avg_value


class BhattaCoefficient(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Bhattacharya coefficient between
    two property vectors

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = "bhatta"

    def __init__(self, property_1, property_2):
        self.min_value = min(
            min(property_1.property_vector), min(property_2.property_vector)
        )
        self.max_value = max(
            max(property_1.property_vector), max(property_2.property_vector)
        )
        super().__init__(property_1, property_2, clean=True)

    def calculate_dissimilarity(self):
        """
        Method that returns the Bhatta coefficient from distance between two vectors
        """
        property_1_discretized = self.property_1.discretize_vector(
            min_value=self.min_value, max_value=self.max_value
        )
        property_2_discretized = self.property_2.discretize_vector(
            min_value=self.min_value, max_value=self.max_value
        )

        return dictances.bhattacharyya(property_1_discretized, property_2_discretized)


class KLDivergence(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Kullback-Leibler divergence between
    two property vectors

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = "kl"

    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_dissimilarity(self):
        """
        Method that returns the KL divergence distance between two vectors
        """
        return dictances.kullback_leibler(
            self.property_1.property_vector_discretized,
            self.property_2.property_vector_discretized,
        )


class PearsonCoefficient(Dissimilarity):
    """
    A Subclass of the Dissimilarity class that calculates the Pearson coefficient between
    two property vectors

    Attributes
    ----------
    property_1 : ProteinProperty object
        Refers to the calculated property of the full protein trajectory
    property_1 : ProteinProperty object
        Refers to the calculated property of the sample protein trajectory
    """

    display_name = "pearson"

    def __init__(self, property_1, property_2):
        super().__init__(property_1, property_2, clean=True)

    def calculate_dissimilarity(self):
        """
        Method that returns the pearson coefficient between two vectors
        """
        return dictances.pearson(
            self.property_1.property_vector_discretized,
            self.property_2.property_vector_discretized,
        )
