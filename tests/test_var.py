from unittest import TestCase
import numpy as np
import nucleoatac.NucleosomeCalling as Nuc
import pyatac.VMat as V
from pyatac.chunkmat2d import BiasMat2D
from pyatac.chunk import ChunkList
from pyatac.bias import InsertionBiasTrack


class Test_variance(TestCase):
    """class for testing variance calculation on background signal

    """
    def setUp(self):
        """ set up class for testing variance calculation for background signal

        """
        bed_list = ChunkList.read('example/example.bed')
        chunk = bed_list[0]
        vmat = V.VMat.open('example/example.VMat')
        biastrack = InsertionBiasTrack(chunk.chrom, chunk.start, chunk.end)
        biastrack.read_track('example/example.Scores.bedgraph.gz')
        biasmat = BiasMat2D(chunk.chrom,chunk.start+200,chunk.end-200,100,250)
        biasmat.makeBiasMat(biastrack)
        self.signaldist = Nuc.SignalDistribution(chunk.start+300,vmat,biasmat,35)
    def test_sd1(self):
        """Make sure variance calculation is close to what is obtained by simulation

        """
        self.signaldist.simulateDist(5000)
        sd1 = np.std(self.signaldist.scores)
        sd2 = self.signaldist.analStd()
        self.assertTrue(abs(sd1-sd2)<0.05*sd1)
    def test_sd2(self):
        """Make sure variance calculation is same as would be obtained through alternate calculation

        """
        var_term = np.sum(self.signaldist.prob_mat*(1-self.signaldist.prob_mat)*self.signaldist.vmat.mat**2)
        tmp = self.signaldist.prob_mat *self.signaldist.vmat.mat
        cov_term = np.sum(np.outer(tmp,tmp))-np.sum(tmp**2)
        sd1 = np.sqrt(self.signaldist.reads * (var_term - cov_term))
        sd2 = self.signaldist.analStd()
        self.assertTrue(abs(sd1-sd2)<0.001*sd1)


