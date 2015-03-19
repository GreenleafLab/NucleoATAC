from unittest import TestCase

import numpy as np
from pyatac.utils import call_peaks
from pyatac.chunkmat2d import FragmentMat2D, BiasMat2D
from pyatac.chunk import ChunkList
from pyatac.bias import InsertionBiasTrack
from pyatac.fragmentsizes import FragmentSizes
from pyatac.tracks import InsertionTrack, Track

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

class Test_ReadMat2D(TestCase):
    """test methods for ReadMat2D"""
    def test_get(self):
        """test get function for fragmentmat"""
        x = FragmentMat2D('chr1',500,1000,0,200)
        x.mat[100,5] = 1
        self.assertTrue(np.array_equal(x.get(start=505,end=507,i_lower=100,i_upper=102),np.array([[1,0],[0,0]])))


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
        isizes = FragmentSizes(i_lower=100,i_upper=200, vals = np.array(range(100,200)))
        self.biasmat.normByInsertDist(isizes)
        a1 = self.biastrack.get(pos = self.biasmat.start -50)
        a2 = self.biastrack.get(pos = self.biasmat.start + 50)
        correct = np.exp(a1+a2)*isizes.get(size = 101)
        self.assertTrue(abs(correct - self.biasmat.mat[1,0])<0.01*correct)

class Test_Ins(TestCase):
    """test that conversion of readmat to insertion gives same result as insertion track"""
    def setUp(self):
        """setup Test_Ins class by making a fragmentlist"""
        bed_list = ChunkList.read('example/example.bed')
        self.chunk = bed_list[0]
    def test_ins_methods(self):
        """test that two methods for getting insertion track give same result"""
        ins1 = InsertionTrack(self.chunk.chrom, self.chunk.start, self.chunk.end)
        ins1.calculateInsertions('example/single_read.bam')
        mat = FragmentMat2D(self.chunk.chrom,self.chunk.start,self.chunk.end,0,100)
        mat.makeFragmentMat(self.fragmentlist)
        ins2 = mat.getIns()
        self.assertTrue(np.array_equal(ins1.get(self.chunk.start+100,self.chunk.start+300),ins2.get(self.chunk.start+100,self.chunk.start+300)))

class Test_Track(TestCase):
    """Test out the Track class in PyATAC"""
    def setUp(self):
        """setup Test_Track class"""
        bed_list = ChunkList.read('example/example.bed')
        self.chunk = bed_list[0]
    def test_read_and_get(self):
        """test the read and get functionality of track class"""
        track = Track(self.chunk.chrom,self.chunk.start,self.chunk.end)
        track.read_track('example/example.Scores.bedgraph.gz')
        val = 1.35994655714
        self.assertTrue(abs(val - track.get(pos = 706661))<0.001)





