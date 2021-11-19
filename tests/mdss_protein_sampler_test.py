import random
import mdss_protein_sampler


def test_calculate_distance():
    class FakeProperty:
        def __init__(self):
            self.avg_value = random.random()

    property_1 = FakeProperty()
    property_2 = FakeProperty()
    distance = mdss_protein_sampler.Distance(property_1, property_2)

    assert distance.distance == property_1.avg_value - property_2.avg_value
