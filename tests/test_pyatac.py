from unittest import TestCase

import nucleoatac.PyATAC as PA
import numpy as np
import bx.bbi.bigwig_file
import pysam

class Test_call_peaks(TestCase):
    """test the call_peaks function in pyatac"""
    def test_call_peaks_case1(self):
        """test case1 for call_peaks"""
        peaks = PA.call_peaks(np.array([1,2,3,2,1,4,1,2,1,0,0]),min_signal=1,sep=3)
        self.assertTrue(np.array_equal(peaks,np.array([2,5])))
    def test_call_peaks_case2(self):
        """test case2 for call_peaks"""
        peaks = PA.call_peaks(np.array([1,2,3,2,1,4,1,2,1,0,0]),min_signal=1,sep=1)
        self.assertTrue(np.array_equal(peaks,np.array([2,5,7])))
    def test_call_peaks_case3(self):
        """test case3 for call_peaks"""
        peaks = PA.call_peaks(np.array([1,2,3,2,1,4,1,2,1,0,0]),min_signal=3,sep=2)
        self.assertTrue(np.array_equal(peaks,np.array([2,5])))

class Test_ReadMat2D(TestCase):
    """test methods for ReadMat2D"""
    def test_updatedMat(self):
        """test updating of readmat"""
        x = PA.ReadMat2D('chrI', 500, 1000, 0, 200)
        read = PA.Read('chrI', 550, 100)
        x.updateMat(read)
        y = np.zeros((200,500))
        y[100,99]=1
        self.assertTrue(np.array_equal(x.get(),y)) 
    def test_get(self):
        """test get function for readmat"""
        x = PA.ReadMat2D('chr1',500,1000,0,200)
        x.mat[100,5] = 1
        self.assertTrue(np.array_equal(x.get(start=505,end=507,i_lower=100,i_upper=102),np.array([[1,0],[0,0]])))

class Test_readlist(TestCase):
    """tests readlist class"""
    def setUp(self):
        """read in readlist"""
        bam = pysam.Samfile('example/single_read.bam','rb')
        self.readlist = PA.ReadList("chrII",706800,707000)
        self.readlist.extract_reads(bam)
    def test_chrom(self):
        """test chromosome of read"""
        self.assertTrue(self.readlist.reads[0].chrom == "chrII")
    def test_left(self):
        """test left position of read"""
        self.assertTrue(self.readlist.reads[0].left ==  706858)
    def test_right(self):
        """test right position of read"""
        self.assertTrue(self.readlist.reads[0].right == 706929)
    def test_insert(self):
        """test insert size of read"""
        self.assertTrue(self.readlist.reads[0].insert == 71)

class Test_BiasMat(TestCase):
    """test construction of BiasMat"""
    def setUp(self):
        """setup Test_BiasMat class with construction of a biasmat"""
        bed_list = PA.read_bed('example/example.bed')
        self.chunk = bed_list[0]
        bias_bwh=bx.bbi.bigwig_file.BigWigFile(open('example/example.Scores.bw','rb'))
        self.biastrack = PA.InsertionBiasTrack(self.chunk.chrom, self.chunk.start, self.chunk.end)
        self.biastrack.read_track(bias_bwh)
        self.biasmat = PA.BiasMat2D(self.chunk.chrom,self.chunk.start+100,self.chunk.end-100,100,200)
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
        isizes = PA.InsertSizes(np.array(range(100,200)),i_lower=100,i_upper=200)
        self.biasmat.normByInsertDist(isizes)
        a1 = self.biastrack.get(pos = self.biasmat.start -50)
        a2 = self.biastrack.get(pos = self.biasmat.start + 50)
        correct = np.exp(a1+a2)*isizes.get(size = 101)
        self.assertTrue(abs(correct - self.biasmat.mat[1,0])<0.01*correct)

class Test_Ins(TestCase):
    """test that conversion of readmat to insertion gives same result as insertion track"""
    def setUp(self):
        """setup Test_Ins class by making a readlist"""
        bed_list = PA.read_bed('example/example.bed')
        self.chunk = bed_list[0]
        self.readlist = PA.ReadList(self.chunk.chrom, self.chunk.start, self.chunk.end)
        self.readlist.addRead(PA.Read(self.chunk.chrom, self.chunk.start + 200, 93))
        self.readlist.addRead(PA.Read(self.chunk.chrom, self.chunk.start + 155, 43))
        self.readlist.addRead(PA.Read(self.chunk.chrom, self.chunk.start + 138, 90))
        self.readlist.addRead(PA.Read(self.chunk.chrom, self.chunk.start + 155, 43))
        self.readlist.addRead(PA.Read(self.chunk.chrom, self.chunk.start + 147, 73))
    def test_ins_methods(self):
        """test that two methods for getting insertion track give same result"""
        ins1 = PA.InsertionTrack(self.chunk.chrom, self.chunk.start, self.chunk.end)
        ins1.calculateInsertions(self.readlist)
        mat = PA.ReadMat2D(self.chunk.chrom,self.chunk.start,self.chunk.end,0,100)
        mat.makeReadMat(self.readlist)
        ins2 = mat.getIns()
        self.assertTrue(np.array_equal(ins1.get(self.chunk.start+100,self.chunk.start+300),ins2.get(self.chunk.start+100,self.chunk.start+300)))

class Test_Track(TestCase):
    """Test out the Track class in PyATAC"""
    def setUp(self):
        """setup Test_Track class"""
        bed_list = PA.read_bed('example/example.bed')
        self.chunk = bed_list[0]
    def test_read_and_get(self):
        """test the read and get functionality of track class"""
        track = PA.Track(self.chunk.chrom,self.chunk.start,self.chunk.end)
        bias_bwh = bx.bbi.bigwig_file.BigWigFile(open('example/example.Scores.bw','rb')) 
        track.read_track(bias_bwh)
        val = 1.35994655714
        self.assertTrue(abs(val - track.get(pos = 706661))<0.001)





