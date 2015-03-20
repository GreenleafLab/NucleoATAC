from unittest import TestCase

import numpy as np
from pyatac.chunkmat2d import FragmentMat2D, BiasMat2D
from pyatac.chunk import ChunkList
from pyatac.bias import InsertionBiasTrack
from pyatac.fragmentsizes import FragmentSizes
from pyatac.tracks import InsertionTrack, Track


class Test_FragmentMat2D(TestCase):
    """test methods for ReadMat2D"""
    def test_get(self):
        """test get function for fragmentmat"""
        x = FragmentMat2D('chr1',500,1000,0,200)
        x.mat[100,5] = 1
        self.assertTrue(np.array_equal(x.get(start=505,end=507,lower=100,upper=102),np.array([[1,0],[0,0]])))


class Test_BiasMat(TestCase):
    """test construction of BiasMat"""
    def setUp(self):
        """setup Test_BiasMat class with construction of a biasmat"""
        bed_list = ChunkList.read('example/example.bed')
        self.chunk = bed_list[0]
        self.biastrack = InsertionBiasTrack(self.chunk.chrom, self.chunk.start, self.chunk.end)
        self.biastrack.read_track('example/example.Scores.bedgraph.gz')
        self.biasmat = BiasMat2D(self.chunk.chrom,self.chunk.start+100,self.chunk.end-100,100,200)
        self.biasmat.makeBiasMat(self.biastrack)
    def test_biasmat1(self):
        """test case1 for biasmat"""
        a1 = self.biastrack.get(pos = self.biasmat.start - 49)
        a2 = self.biastrack.get(pos = self.biasmat.start + 50)
        correct = np.exp(a1+a2)
        self.assertTrue(abs(correct - self.biasmat.mat[0,0])<0.01*correct)
    def test_biasmat2(self):
        """test case2 for biasmat"""
        a1 = self.biastrack.get(pos = self.biasmat.start + 145)
        a2 = self.biastrack.get(pos = self.biasmat.start + 295)
        correct = np.exp(a1+a2)
        self.assertTrue(abs(correct - self.biasmat.mat[51,220]) < 0.01*correct)
    def test_normByInsertDist(self):
        """test that normalization by insert distribution works as expected"""
        isizes = FragmentSizes(lower=100,upper=200, vals = np.array(range(100,200)))
        self.biasmat.normByInsertDist(isizes)
        a1 = self.biastrack.get(pos = self.biasmat.start -50)
        a2 = self.biastrack.get(pos = self.biasmat.start + 50)
        correct = np.exp(a1+a2)*isizes.get(size = 101)
        self.assertTrue(abs(correct - self.biasmat.mat[1,0])<0.01*correct)


