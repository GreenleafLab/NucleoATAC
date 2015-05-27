"""
General tools for dealing with ATAC-Seq data using Python.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import numpy as np
import matplotlib.pyplot as plt
from pyatac.bedgraph import BedGraphFile
from pyatac.chunk import Chunk
from pyatac.utils import smooth
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
from fragments import getInsertions, getStrandedInsertions
from pyatac.seq import get_sequence, seq_to_mat, complement

class Track(Chunk):
    """Generic class for various types of signal tracks"""
    def __init__(self, chrom, start, end, name = "track", vals=None , log = False):
        Chunk.__init__(self, chrom, start, end, name = name)
        self.log = log
        if vals is None:
            self.vals = None
        elif len(vals) == self.length():
            self.vals = vals
        else:
            raise Exception("Input vals must be of length as set by start and end!")
    def assign_track(self, vals, start = None, end = None):
        """Assign values to track"""
        if start:
            self.start = start
        if end:
            self.end = end
        if len(vals)!= self.end - self.start:
            raise Exception("The values being assigned to track do not \
                    span the start to end of the track")
        self.vals = vals
    def write_track(self, handle, start = None, end = None, vals = None, write_zero = True):
        """Write track to output file handle

        If vals are specified use those values
        Othersise use self.vals
        """
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        if vals is None:
            vals=self.vals
        if len(vals)!=self.end-self.start:
            print len(vals),self.end-self.start
            raise Exception("Error! Inconsistency between length of \
            values and start/end values")
        prev_value = None
        start_range = 0
        output = ""
        for i in range(len(vals)):
            if vals[i] == prev_value:
                pass
            elif np.isnan(vals[i]):
                prev_value = vals[i]
            elif prev_value is not None and not np.isnan(prev_value):
                if write_zero or prev_value!=0:
                    output += "\t".join(map(str,[self.chrom,start_range,start+i,prev_value]))+"\n"
                start_range = start + i
                prev_value = vals[i]
            else:
                start_range = start + i
                prev_value = vals[i]
        if prev_value==0:
            if write_zero:
                output += "\t".join(map(str,[self.chrom,start_range,end,prev_value]))+"\n"
        elif not np.isnan(prev_value):
            output += "\t".join(map(str,[self.chrom,start_range,end,prev_value]))+"\n"
        handle.write(output)
    def read_track(self, bedgraph, start = None, end = None, empty = np.nan, flank = None):
        """Read track values from BigWig file handle"""
        if start:
            self.start = start
        if end:
            self.end = end
        if flank:
            self.start = self.start - flank
            self.end = self.end + flank
        handle = BedGraphFile(bedgraph)
        self.vals = handle.read(self.chrom, self.start,
                                     self.end, empty = empty)
        handle.close()
    def log(self, pseudo = 1):
        """Log values.  Add psuedo count so values don't equal 0 before logging"""
        if self.log:
            print "Logging a track that is already log..."
        adjusted = self.vals + pseudo
        self.vals = np.log(adjusted)
        self.log = True
    def exp(self):
        """Take exponent of values"""
        if not self.log:
            print "taking exponent of a non-logged track..."
        self.vals = np.exp(self.vals)
        self.log = False
    def smooth_track(self,  window_len, window='flat', sd = None,
                     mode = 'valid', norm = True):
        """ smoothing of track"""
        self.smoothed=True
        self.vals = smooth(self.vals, window_len, window = window, sd = sd,
                           mode = mode, norm = norm)
        if mode == 'valid':
            self.start = self.start + window_len/2
            self.end = self.end - window_len/2
    def get(self, start = None, end = None, pos = None):
        """Obtain value of track at particular interval or position"""
        if pos:
            try:
                return self.vals[pos-self.start]
            except:
                raise Exception("Looks like position given doesn't match track")
        else:
            if start is None:
                start = self.start
            if end is None:
                end = self.end
            x1 = start - self.start
            x2 = end - self.start
            try:
                return self.vals[x1:x2]
            except:
                raise Exception("Looks like dimensions from get probaby don't match track, or there are no vals in track")
    def plot(self, name = None, start = None, end = None, filename=None):
        """plot the values """
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        if name is None:
            name = self.name
        values = self.get(start,end)
        fig = plt.figure()
        plt.plot(range(start,end),values)
        plt.xlabel(self.chrom)
        plt.ylabel(name)
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            #Also save text output!
            filename2 = ".".join(filename.split(".")[:-1]+['txt'])
            out = np.vstack((range(start,end),values))
            np.savetxt(filename2,out,delimiter="\t")
        else:
            fig.show()
    def slop(self, chromDict, up = 0, down = 0, new = False):
        """Modification of slop method for Chunk class to check if vals are set"""
        if self.vals is None:
            out = Chunk.slop(self, chromDict, up = up, down = down, new = new)
            return out
        else:
            raise Exception("Cannot slop Track if vals are set")



class InsertionTrack(Track):
    """Class for getting and storing insertion positions"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "insertions")
    def calculateInsertions(self, bamfile, flank = 0, lower = 0, upper = 2000, atac = True):
        """Compute inserion track"""
        self.start = self.start - flank
        self.end = self.end + flank
        self.vals = getInsertions(bamfile, self.chrom, self.start, self.end, lower, upper, atac)
    def calculateStrandedInsertions(self, bamfile, flank =0, lower = 0, upper = 2000, atac = True):
        """Compute inserion track for plus and minus strands separately"""
        self.start = self.start - flank
        self.end = self.end + flank
        self.plus, self.minus = getStrandedInsertions(bamfile, self.chrom, self.start, self.end, lower, upper, atac)
        self.vals = self.plus + self.minus
    def getInsertionSequences(self, fasta, nucleotides = ["C","G","A","T"], up = 10, down = 10):
        """Get sequence content at insertions"""
        mat = np.zeros((len(nucleotides), up + down +1))
        if np.sum(self.vals) == 0:
            return mat
        offset = max(up,down)
        seq_chunk = Chunk(self.chrom, self.start - offset, self.end + offset)
        sequence = get_sequence(seq_chunk, fasta)
        seq_mat = seq_to_mat(sequence, nucleotides)
        for i in range(self.length()):
            mat += self.vals[i] * seq_mat[:,(offset + i - up):(offset + i + down + 1)]
        return mat
    def getStrandedInsertionSequences(self, fasta, nucleotides = ["C","G","A","T"], up = 10, down = 10):
        """Get sequence content at insertions, taking into account strand"""
        mat = np.zeros((len(nucleotides), up + down +1))
        if np.sum(self.vals) == 0:
            return mat
        offset = max(up,down)
        seq_chunk = Chunk(self.chrom, self.start - offset, self.end + offset)
        sequence = get_sequence(seq_chunk, fasta)
        minus_sequence = complement(sequence)
        seq_mat = seq_to_mat(sequence, nucleotides)
        minus_seq_mat = seq_to_mat(minus_sequence, nucleotides)
        for i in range(self.length()):
            mat += self.plus[i] * seq_mat[:,(offset + i - up):(offset + i + down + 1)]
            mat += self.minus[i] * np.fliplr(minus_seq_mat[:,(offset + i - down):(offset + i + up + 1)])
        return mat



class CoverageTrack(Track):
    """Class for computing read center converage"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "coverage")
    def calculateCoverage(self, mat, lower, upper, window_len):
        """Compute coverage of fragment centers using flat window"""
        offset=self.start-mat.start-(window_len/2)
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate coverage with desired window")
        lower=lower-mat.lower
        upper=upper-mat.lower
        if offset!=0:
            collapsed = np.sum(mat.mat[lower:upper,offset:-offset],axis=0)
        else:
            collapsed = np.sum(mat.mat[lower:upper,],axis=0)
        self.vals = smooth(collapsed, window_len, window="flat",
                            mode='valid',norm=False)
    def calculateCoverageSmooth(self,mat,lower,upper,window_len,sd):
        """Compute coverage of fragment centers using gaussia window"""
        offset=self.start-mat.start-(window_len/2)
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate coverage with desired window")
        lower=lower-mat.lower
        upper=upper-mat.lower
        if offset!=0:
            collapsed = np.sum(mat.mat[lower:upper,offset:-offset],axis=0)
        else:
            collapsed = np.sum(mat.mat[lower:upper,],axis=0)
        self.vals = smooth(collapsed, window_len, sd= sd, window="gaussian",
                            mode='valid',norm=False)






