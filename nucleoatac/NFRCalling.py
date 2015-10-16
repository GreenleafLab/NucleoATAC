"""
Script with classes and functions for nucleosome calling.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import numpy as np
import pysam
from pyatac.tracks import Track, InsertionTrack
from pyatac.chunk import Chunk
from pyatac.utils import read_chrom_sizes_from_fasta
from pyatac.bias import PWM, InsertionBiasTrack
from copy import copy


class NFR(Chunk):
    """Class for storing information about a single NFR region"""
    def __init__(self, left, right, nfrtrack):
        self.chrom = nfrtrack.chrom
        self.start = left
        self.end = right
        self.strand = "*"
        self.occ = np.mean(nfrtrack.occ.get(left,right))
        self.min_upper = np.min(nfrtrack.occ_upper.get(left,right))
        self.ins_density = np.mean(nfrtrack.ins.get(left,right))
        self.bias_density = np.mean(nfrtrack.bias.get(left, right, log = False))
    def asBed(self):
        out = "\t".join(map(str,[self.chrom, self.start, self.end, self.occ, self.min_upper, self.ins_density, self.bias_density]))
        return out 
    def write(self, handle):
        handle.write(self.asBed() + "\n")


class NFRParameters:
    """Class for storing parameters related to nucleosome calling"""
    def __init__(self, occ_track,calls, ins_track = None, bam = None, max_occ = 0.25, max_occ_upper = 0.25, fasta = None, pwm = None):
        self.bam = bam
        self.ins_track = ins_track
        self.occ_track = occ_track
        self.calls = calls
        self.max_occ = max_occ
        self.max_occ_upper = max_occ_upper
        self.fasta = fasta
        if fasta is not None:
            self.pwm = PWM.open(pwm)
            self.chrs = read_chrom_sizes_from_fasta(fasta)
  


class NFRChunk(Chunk):
    """Class for storing and determining collection of nfr positions
    """
    def __init__(self, chunk):
        self.start = chunk.start
        self.end = chunk.end
        self.chrom = chunk.chrom
        self.nfrs = []
    def initialize(self, parameters):
        self.params = parameters
    def getOcc(self):
        """gets occupancy track-- reads in from bedgraph"""
        self.occ = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ.read_track(self.params.occ_track)
        #lower_file = self.params.occ_track[:-11] + 'lower_bound.bedgraph.gz'
        #self.occ_lower = Track(self.chrom,self.start,self.end,"Occupancy")
        #self.occ_lower.read_track(lower_file)
        upper_file = self.params.occ_track[:-11] + 'upper_bound.bedgraph.gz'
        self.occ_upper = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ_upper.read_track(upper_file)
    def getIns(self):
        """gets insertion track-- reads in from bedgraph"""
        if self.params.ins_track is None:
            self.ins = InsertionTrack(self.chrom, self.start, self.end)
            self.ins.calculateInsertions(self.params.bam)
        else:
            self.ins = Track(self.chrom,self.start,self.end,"Insertion")
            self.ins.read_track(self.params.ins_track)
    def getBias(self):
        """get bias"""
        self.bias = InsertionBiasTrack(self.chrom, self.start,
                                  self.end, log = True)
        if self.params.fasta is not None:
            self.bias.computeBias(self.params.fasta, self.params.chrs, self.params.pwm)
    def findNFRs(self):
        """find NFR regions"""
        region = np.ones(self.length())
        tbx = pysam.TabixFile(self.params.calls)
        nucs = []
        if self.chrom in tbx.contigs:
            for row in tbx.fetch(self.chrom, self.start, self.end, parser = pysam.asTuple()):
                nucs.append(int(row[1]))
        for j in xrange(1,len(nucs)):
            left = nucs[j-1] + 73
            right = nucs[j] - 72
            if right <= left:
                continue
            candidate = NFR(left, right, self)
            if candidate.min_upper < self.params.max_occ_upper and candidate.occ < self.params.max_occ:
                self.nfrs.append(candidate)
    def process(self, params):
        """wrapper to carry out all methods needed to call nucleosomes and nfrs"""
        self.initialize(params)
        self.getOcc()
        self.getIns()
        self.getBias()
        self.findNFRs()
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self,name)

