from unittest import TestCase

import numpy as np
from pyatac.chunkmat2d import FragmentMat2D
from pyatac.chunk import ChunkList
from pyatac.bias import InsertionBiasTrack
from pyatac.tracks import InsertionTrack, Track


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
        mat.makeFragmentMat('example/single_read.bam')
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





