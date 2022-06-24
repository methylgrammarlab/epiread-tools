import unittest
from em_utils import find_intersection
import numpy as np

class MyTestCase(unittest.TestCase):

    def test_intersection(self):
        #no intersection
        assert not len(find_intersection(np.array([1,3,100,110]), 50, 55))
        #two intersections
        assert find_intersection(np.array([1,3,100,110]), 2, 105) == [(2, 3), (100, 105)]



# if __name__ == '__main__':
#     unittest.main()
