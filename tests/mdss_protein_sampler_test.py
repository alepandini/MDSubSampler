import random
from mdss_protein_sampler import RandomSampler
import mdss_protein_sampler
import pytest


def test_random_sampler_sample_has_expected_length():
    frame_list = list(range(1000))
    random_sampler = mdss_protein_sampler.RandomSampler(frame_list, seed_number=1999)
    sample_list = random_sampler.sample(100)

    assert len(sample_list) == 100


def test_random_sampler_sample_returns_subset():
    frame_list = list(range(1000))
    random_sampler = mdss_protein_sampler.RandomSampler(frame_list, seed_number=1999)
    sample_list = random_sampler.sample(100)

    assert all(x in frame_list for x in sample_list)


def test_calculate_distance():
    class FakeProperty:
        def __init__(self):
            self.avg_value = random.random()

    property_1 = FakeProperty()
    property_2 = FakeProperty()
    distance = mdss_protein_sampler.Distance(property_1, property_2)

    assert distance.distance == property_1.avg_value - property_2.avg_value


@pytest.mark.parametrize(
    "mocked_function_name, distance_subclass",
    [
        ("dictances.bhattacharyya", mdss_protein_sampler.BhattaDistance),
        ("dictances.kullback_leibler", mdss_protein_sampler.KLDiverDistance),
        ("dictances.pearson", mdss_protein_sampler.PearsonDictDistance),
    ],
)
def test_bhattacharyya_distance(mocker, mocked_function_name, distance_subclass):
    class FakeProperty:
        def __init__(self):
            self.property_vector_discretized = random.random()

    mocked_function = mocker.patch(mocked_function_name)

    property_1 = FakeProperty()
    property_2 = FakeProperty()
    distance_subclass(property_1, property_2)

    mocked_function.assert_called_once_with(
        property_1.property_vector_discretized, property_2.property_vector_discretized
    )
