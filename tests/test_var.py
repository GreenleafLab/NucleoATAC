from unittest import TestCase

import nucleoatac.NucleosomeCalling as Nuc
import nucleoatac.VMat as V
import nucleoatac.PyATAC as PA
import numpy as np
import bx.bbi.bigwig_file

class Test_variance(TestCase):
    """class for testing variance calculation on background signal

    """
    def setUp(self):
        """ set up class for testing variance calculation for background signal

        """
        bed_list = PA.read_bed('example/example.bed')
        chunk = bed_list[0]
        vmat = V.VMat.open('example/example.VMat')
        bias_bwh=bx.bbi.bigwig_file.BigWigFile(open('example/example.Scores.bw','rb'))
        biastrack = PA.InsertionBiasTrack(chunk.chrom, chunk.start, chunk.end)
        biastrack.read_track(bias_bwh)
        biasmat = PA.BiasMat2D(chunk.chrom,chunk.start+200,chunk.end-200,100,250)
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


