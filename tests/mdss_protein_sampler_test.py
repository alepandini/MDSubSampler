import random
import mdss_protein_sampler

import dictances


def test_calculate_distance():
    class FakeProperty:
        def __init__(self):
            self.avg_value = random.random()

    property_1 = FakeProperty()
    property_2 = FakeProperty()
    distance = mdss_protein_sampler.Distance(property_1, property_2)

    assert distance.distance == property_1.avg_value - property_2.avg_value


def test_bhattacharyya_distance(mocker):
    class FakeProperty:
        def __init__(self):
            self.property_vector_discretized = random.random()

    mocker.patch("dictances.bhattacharyya")

    property_1 = FakeProperty()
    property_2 = FakeProperty()
    distance = mdss_protein_sampler.BhattaDistance(property_1, property_2)

    dictances.bhattacharyya.assert_called_once_with(
        property_1.property_vector_discretized, property_2.property_vector_discretized
    )
