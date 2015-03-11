from unittest import TestCase

import nucleoatac.PyATAC as PA
import nucleoatac.NucleosomeCalling as Nuc
import nucleoatac.VMat as V
import numpy as np
import pysam

class Test_xcor(TestCase):
    """class to test occupancy"""
    def setUp(self):
        """setup Test_occupancy class by establishing parameters"""
        bed_list = PA.read_bed('example/example.bed')
        self.chunk = bed_list[0]
        bam = pysam.Samfile('example/example.bam','rb')
        self.vmat = V.VMat.open('example/example.VMat')
        readlist = PA.ReadList(self.chunk.chrom, self.chunk.start - self.vmat.w, self.chunk.end + self.vmat.w) 
        readlist.extract_reads(bam)
        self.vmat = V.VMat.open('example/example.VMat')
        self.mat = PA.ReadMat2D(self.chunk.chrom,self.chunk.start-self.vmat.w,self.chunk.end+self.vmat.w,self.vmat.i_lower,self.vmat.i_upper) 
        self.mat.makeReadMat(readlist)
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

