import pytest
from epiread_tools.em_utils import find_intersection
import numpy as np


def test_intersection(self):
    #no intersection
    assert not len(find_intersection(np.array([1,3,100,110]), 50, 55))
    #two intersections
    assert find_intersection(np.array([1,3,100,110]), 2, 105) == [(2, 3), (100, 105)]

