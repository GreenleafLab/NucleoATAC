from unittest import TestCase

import nucleoatac.PyATAC as PA
import nucleoatac.Occupancy as Occ
import numpy as np

class Test_occupancy(TestCase):
    """class to test occupancy"""
    def setUp(self):
        """setup Test_occupancy class by establishing parameters"""
        self.insert_dist = Occ.InsertDistribution(0,3)
        self.insert_dist.nfr_fit = PA.InsertSizes(np.array([0.5, 0.49, 0.01 ]),0,3)
        self.insert_dist.nuc_fit = PA.InsertSizes(np.array([0.01, 0.49,0.5]),0,3) 
        self.params = Occ.OccupancyCalcParams(0,3,self.insert_dist)
    def test_occupancy_calc1(self):
        """test occupancy class for a zero occupancy situation"""
        occ = Occ.calculateOccupancy(np.array([1,0,0]),self.params)
        self.assertTrue(occ == 0)
    def test_occupancy_calc2(self):
        """test occupancy for a 50% occupancy situation"""
        occ = Occ.calculateOccupancy(np.array([1,1,1]),self.params)         
        self.assertTrue(occ == 0.5)
    def test_occupancy_calc3(self):
        """test average occupancy for when you have 40 random reads,
             which should be 75% from nuc dist
        """
        results = np.zeros(100)
        for i in range(100):
            a = np.random.multinomial(10,self.insert_dist.nfr_fit.get())
            b = np.random.multinomial(30,self.insert_dist.nuc_fit.get())
            results[i] = Occ.calculateOccupancy(a+b,self.params)
        self.assertTrue(abs(np.mean(results)-0.75)<0.1)

