from unittest import TestCase
import numpy as np

import nucleoatac.NucleosomeCalling as Nuc
import pyatac.VMat as V
from pyatac.chunkmat2d import FragmentMat2D
from pyatac.chunk import ChunkList



class Test_xcor(TestCase):
    """class to test occupancy"""
    def setUp(self):
        """setup Test_occupancy class by establishing parameters"""
        bed_list = ChunkList.read('example/example.bed')
        self.chunk = bed_list[0]
        self.vmat = V.VMat.open('example/example.VMat')
        self.vmat = V.VMat.open('example/example.VMat')
        self.mat = FragmentMat2D(self.chunk.chrom,self.chunk.start-self.vmat.w,self.chunk.end+self.vmat.w,self.vmat.lower,self.vmat.upper)
        self.mat.makeFragmentMat('example/example.bam')
        self.signal = Nuc.SignalTrack(self.chunk.chrom,self.chunk.start,self.chunk.end)
        self.signal.calculateSignal(self.mat, self.vmat)
    def test_signal_calc1(self):
        """test signal calc"""
        a=np.sum(self.mat.get(start = self.chunk.start, end = self.chunk.start + self.vmat.w*2 + 1)*self.vmat.mat)
        b=self.signal.get(pos=self.chunk.start + self.vmat.w)
        self.assertTrue(abs(a-b)<0.0001)
    def test_signal_calc2(self):
        """test signal calc"""
        a=np.sum(self.mat.get(start = self.chunk.start+100, end = self.chunk.start+ 100 + self.vmat.w*2 + 1)*self.vmat.mat)
        b=self.signal.get(pos=self.chunk.start+ 100 + self.vmat.w)
        self.assertTrue(abs(a-b)<0.0001)

