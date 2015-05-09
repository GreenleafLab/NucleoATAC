"""
Script with classes and functions for nucleosome calling.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import numpy as np
import pysam
from pyatac.tracks import Track
from pyatac.chunk import Chunk
from copy import copy

class NFR(Chunk):
    """Class for storing information about a single NFR region"""
    def __init__(self, left, right, nfrtrack):
        self.chrom = nfrtrack.chrom
        self.start = left
        self.end = right
        self.strand = "*"
        self.occ = np.mean(nfrtrack.occ.vals[(left - nfrtrack.occ.start):(right -nfrtrack.occ.start)])
        self.min_upper = np.min(nfrtrack.occ_upper.vals[(left - nfrtrack.occ.start):(right -nfrtrack.occ.start)])
        self.min_lower = np.min(nfrtrack.occ_lower.vals[(left - nfrtrack.occ.start):(right -nfrtrack.occ.start)])
    def asBed(self):
        out = "\t".join(map(str,[self.chrom, self.start, self.end, self.occ, self.min_lower, self.min_upper]))
        return out 
    def write(self, handle):
        handle.write(self.asBed() + "\n")



class NFRParameters:
    """Class for storing parameters related to nucleosome calling"""
    def __init__(self, occ_track, calls, max_occ = 0.1, max_occ_upper = 0.25, max_nfr_gap = 10, min_nfr_len = 10):
        self.occ_track = occ_track
        self.calls = calls
        self.max_nfr_gap = max_nfr_gap
        self.max_occ = max_occ
        self.max_occ_upper = max_occ_upper
        self.min_nfr_len = min_nfr_len

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
        """gets occupancy track-- either reads in from bw handle given, or makes new"""
        self.occ = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ.read_track(self.params.occ_track)
        lower_file = self.params.occ_track[:-11] + 'lower_bound.bedgraph.gz'
        self.occ_lower = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ_lower.read_track(lower_file)
        upper_file = self.params.occ_track[:-11] + 'upper_bound.bedgraph.gz'
        self.occ_upper = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ_upper.read_track(upper_file)
    def findNFRs(self):
        """find NFR regions"""
        region = np.ones(self.length())
        tbx = pysam.TabixFile(self.params.calls)
        for row in tbx.fetch(self.chrom, self.start, self.end, parser = pysam.asTuple()):
            region[max(int(row[1]) - self.start - 73,0):min(int(row[2])-self.start + 73,self.end-self.start)] = 0
        tbx.close()
        tmp = copy(self.occ.vals)
        tmp[np.where(np.isnan(tmp))] = 1
        tmp2 = copy(self.occ_upper.vals)
        tmp2[np.where(np.isnan(tmp2))] = 1
        region[np.where(tmp > self.params.max_occ)] = 0
        region[np.where(tmp2 > self.params.max_occ_upper)] = 0
        starts = [i for i in xrange(1,len(region)-1) if (region[i] == 1 and region[i-1] == 0)]
        if len(starts) == 0:
            return
        ends = [i for i in xrange(starts[0],len(region)) if (region[i] ==0 and region[i-1] ==1)]
        if len(ends) == 0:
            return
        prev_end = ends[0]
        prev_start = starts[0]
        for i in xrange(1,len(starts)+1):
            if i >= len(ends):
                if (prev_end - prev_start) > self.params.min_nfr_len:
                    self.nfrs.append(NFR(self.start + prev_start, self.start + prev_end, self))
                return
            if starts[i] < prev_end + self.params.max_nfr_gap:
                prev_end = ends[i]
            else:
                if (prev_end - prev_start) > self.params.min_nfr_len:
                    self.nfrs.append(NFR(self.start + prev_start, self.start + prev_end, self))
                prev_start = starts[i]
                prev_end = ends[i]
    def process(self, params):
        """wrapper to carry out all methods needed to call nucleosomes and nfrs"""
        self.initialize(params)
        self.getOcc()
        self.findNFRs()
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self,name)

