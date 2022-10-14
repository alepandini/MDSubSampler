import mdss_dissimilarity
import random
import pytest


def test_calculate_distance():
    class FakeProperty:
        def __init__(self):
            self.avg_value = random.random()

    property_1 = FakeProperty()
    property_2 = FakeProperty()
    distance = mdss_dissimilarity.Dissimilarity(property_1, property_2)

    assert distance.distance == property_1.avg_value - property_2.avg_value


@pytest.mark.parametrize(
    "mocked_function_name, distance_subclass",
    [
        ("dictances.bhattacharyya", mdss_dissimilarity.Bhattacharya),
        ("dictances.kullback_leibler", mdss_dissimilarity.KullbackLeibler),
        ("dictances.pearson", mdss_dissimilarity.Pearson),
    ],
)
def test_bhattacharyya_distance(mocker, mocked_function_name, distance_subclass):
    class FakeProperty:
        def __init__(self):
            self.property_distribution_dict = random.random()

    mocked_function = mocker.patch(mocked_function_name)

    property_1 = FakeProperty()
    property_2 = FakeProperty()
    distance_subclass(property_1, property_2)

    mocked_function.assert_called_once_with(
        property_1.property_distribution_dict, property_2.property_distribution_dict
    )
