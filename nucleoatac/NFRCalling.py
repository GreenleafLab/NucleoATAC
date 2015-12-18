"""
Script with classes and functions for nucleosome calling.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import numpy as np
import pysam
from pyatac.tracks import Track, InsertionTrack, CoverageTrack
from pyatac.chunk import Chunk, _chunkCompare
from pyatac.utils import read_chrom_sizes_from_fasta, call_peaks
from pyatac.bias import PWM, InsertionBiasTrack
from pyatac.chunkmat2d import FragmentMat2D, BiasMat2D
from copy import copy, deepcopy


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
        self.nfr_cov = np.mean(nfrtrack.cov.get(left,right))
    def asBed(self):
        out = "\t".join(map(str,[self.chrom, self.start, self.end, self.occ, self.min_upper, self.ins_density, self.bias_density, self.nfr_cov]))
        return out 
    def write(self, handle):
        handle.write(self.asBed() + "\n")


class NFRParameters:
    """Class for storing parameters related to nucleosome calling"""
    def __init__(self, occ_track,calls, ins_track = None, bam = None, max_occ = 0.25, max_occ_upper = 0.5,
                     fasta = None, pwm = None, upper = 251, flank = 60):
        self.bam = bam
        self.ins_track = ins_track
        self.occ_track = occ_track
        self.calls = calls
        self.max_occ = max_occ
        self.max_occ_upper = max_occ_upper
        self.fasta = fasta
        self.upper = upper
        self.flank = flank
        self.window = flank * 2 +1
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
        self.smooth_ins = deepcopy(self.ins)
        self.smooth_ins.smooth_track(window_len = self.params.window, window = "gaussian", sd = self.params.flank/3, mode='same')
    def getBias(self):
        """get bias"""
        self.bias = InsertionBiasTrack(self.chrom, self.start,
                                  self.end, log = True)
        if self.params.fasta is not None:
            self.bias.computeBias(self.params.fasta, self.params.chrs, self.params.pwm)
    def getNFRCov(self):
        """Get read coverage for regions"""
        self.cov = CoverageTrack(self.chrom, self.start, self.end)
        mat = FragmentMat2D(self.chrom, self.start - self.params.flank,
                                 self.end + self.params.flank, 0, self.params.upper)
        mat.makeFragmentMat(self.params.bam)
        self.cov.calculateCoverage(mat, 0, self.params.upper, self.params.window)
        self.cov.smooth_track(window="gaussian",window_len = self.params.window, sd = self.params.flank/3, mode='same') 
        self.cov.vals *= (1 -self.occ.vals )
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
    def findNFRs2(self):
        tbx = pysam.TabixFile(self.params.calls)
        nucs = []
        if self.chrom in tbx.contigs:
            for row in tbx.fetch(self.chrom, self.start, self.end, parser = pysam.asTuple()):
                if self.occ.get(pos = int(row[1])) > self.params.max_occ:
                    nucs.append(int(row[1]))
        nucs = np.array(nucs)
        candidates = call_peaks(self.cov.get(), sep = self.params.flank, min_signal = 1)
        candidates = candidates[np.where(self.occ.vals[candidates] < self.params.max_occ)] 
        #candidates = candidates[np.where(self.occ_upper.vals[candidates] < self.params.max_occ_upper)] 
        candidates += self.start
        previous_bounds = (-1,-1)
        for candidate in candidates:
            left_nuc = max(np.append(nucs[np.where(nucs < candidate )]+73,self.start))
            right_nuc = max(np.append(nucs[np.where(nucs > candidate )]-73,self.end))
            try:
                left_occ_bound = np.where(self.occ.get(end = candidate) > self.params.max_occ)[0][-1] + 1 + self.start
            except IndexError:
                left_occ_bound = self.start
            try:
                right_occ_bound = np.where(self.occ.get(start = candidate) > self.params.max_occ)[0][0] - 1 + candidate
            except IndexError:
                right_occ_bound = self.end
            try:
                left_cov_bound = np.where(self.cov.get(end = candidate) < self.cov.get(pos = candidate)*0.25)[0][-1] + 1 + self.start
            except IndexError:
                left_cov_bound = self.start
            try:
                right_cov_bound = np.where(self.cov.get(start = candidate) < self.cov.get(pos=candidate) *0.25)[0][0] - 1 + candidate
            except IndexError:
                right_cov_bound = self.end
            left = max(left_occ_bound,left_nuc, left_cov_bound)
            right = min(right_occ_bound,right_nuc,right_cov_bound)
            if right > left and (left,right) != previous_bounds:
                previous_bounds = (left,right)
                nfr_candidate = NFR(left,right,self)
                if nfr_candidate.min_upper < self.params.max_occ_upper:
                    self.nfrs.append(nfr_candidate)
        #merge overlapping
        if len(self.nfrs) > 1:
            self.nfrs.sort(cmp = _chunkCompare)
            i = 0
            merged = []
            previous_left = self.nfrs[0].start
            previous_right = self.nfrs[0].end
            for i in range(1,len(self.nfrs)):
                if self.nfrs[i].start <= previous_right:
                    previous_right = max(self.nfrs[i].end, previous_right)
                else:
                    merged.append(NFR(previous_left, previous_right, self))
                    previous_left = self.nfrs[i].start
                    previous_right = self.nfrs[i].end
            merged.append(NFR(previous_left, previous_right,self))
            self.nfrs = merged
    def process(self, params):
        """wrapper to carry out all methods needed to call nucleosomes and nfrs"""
        self.initialize(params)
        self.getOcc()
        self.getIns()
        self.getBias()
        self.getNFRCov()
        self.findNFRs2()
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self,name)

