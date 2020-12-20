from unittest import TestCase
import numpy as np
import matplotlib
matplotlib.use('agg')

from pyatac.utils import read_chrom_sizes_from_fasta
import nucleoatac.NucleosomeCalling as Nuc
import pyatac.VMat as V
from pyatac.chunkmat2d import FragmentMat2D
from pyatac.chunk import ChunkList
from pyatac.fragmentsizes import FragmentSizes
from pyatac.bias import PWM


class Test_NuclesomeCalling(TestCase):
    """class to test nucleosome calling"""
    def setUp(self):
        """setup Test_occupancy class by establishing parameters"""
        bed_list = ChunkList.read('example/example.bed')
        self.chunk = bed_list[0]
        self.vmat = V.VMat.open('example/example_results/example.VMat')
        #chrs = read_chrom_sizes_from_fasta('example/sacCer3.fa')
        #pwm = PWM.open('Human')
        #chunks = ChunkList.read(
        #    'example/example.bed',
        #    chromDict = chrs,
        #    min_offset = self.vmat.mat.shape[1] + self.vmat.upper//2 + max(pwm.up,pwm.down) + 120//2,
        #    min_length = 120 * 2
        #)
        #chunks.slop(chrs, up = 120//2, down = 120//2)
        #chunks.merge()
        #self.chunk = chunks[0]
        fragment_dist = FragmentSizes.open('example/example_results/example.fragmentsizes.txt')
        self.params = Nuc.NucParameters(vmat = self.vmat, fragmentsizes = fragment_dist,
            bam = 'example/example.bam', fasta = 'example/sacCer3.fa', pwm = 'Human',
            occ_track = 'example/example_results/example.occ.bedgraph.gz',
            sd = 25, nonredundant_sep = 120, redundant_sep = 25,
            min_z = 3, min_lr = 0 , atac = True)
        nuc = Nuc.NucChunk(self.chunk)
        nuc.process(self.params)
        self.out = {
            'nucpos' : [nuc.nuc_collection[i] for i in sorted(nuc.nonredundant)],
            'nucpos.redundant' : [nuc.nuc_collection[i] for i in sorted(nuc.redundant)],
            'nucleoatac_signal' : nuc.norm_signal,
            'nucleoatac_raw' : nuc.nuc_signal,
            'nucleoatac_background' : nuc.bias,
            'nucleoatac_signal.smooth' : nuc.smoothed
        }
    def test_calling(self):
        """test nucleosome positions"""
        self.assertTrue(len(self.out['nucpos']) == 3)
        self.assertTrue(self.out['nucpos'][0].start == 706773)
        self.assertTrue(abs(self.out['nucpos'][0].z - 8.2365) < 0.001)
        self.assertTrue(self.out['nucpos'][0].nfr_cov == 6.0)
        self.assertTrue(self.out['nucpos'][0].nuc_cov == 54.0)
        self.assertTrue(abs(self.out['nucpos'][0].nuc_signal - 15.3150) < 0.001)
        self.assertTrue(abs(self.out['nucpos'][0].norm_signal - 6.4621) < 0.001)
        self.assertTrue(abs(self.out['nucpos'][0].smoothed - 3.4974) < 0.001)
        self.assertTrue(abs(self.out['nucpos'][0].lr - 22.3557) < 0.001)
    def test_signal(self):
        """test signal output"""
        self.assertTrue(len(self.out['nucleoatac_signal'].vals) == 1093)
