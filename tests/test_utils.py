from unittest import TestCase

import numpy as np
from pyatac.utils import call_peaks

class Test_call_peaks(TestCase):
    """test the call_peaks function in pyatac"""
    def test_call_peaks_case1(self):
        """test case1 for call_peaks"""
        peaks = call_peaks(np.array([1,2,3,2,1,4,1,2,1,0,0]),min_signal=1,sep=3)
        self.assertTrue(np.array_equal(peaks,np.array([2,5])))
    def test_call_peaks_case2(self):
        """test case2 for call_peaks"""
        peaks = call_peaks(np.array([1,2,3,2,1,4,1,2,1,0,0]),min_signal=1,sep=1)
        self.assertTrue(np.array_equal(peaks,np.array([2,5,7])))
    def test_call_peaks_case3(self):
        """test case3 for call_peaks"""
        peaks = call_peaks(np.array([1,2,3,2,1,4,1,2,1,0,0]),min_signal=3,sep=2)
        self.assertTrue(np.array_equal(peaks,np.array([2,5])))


