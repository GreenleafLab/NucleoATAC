#### Import needed modules #####
import pyatac.seq as seq
import numpy as np
from pysam import AlignmentFile, AlignedSegment
from pysam import FastaFile

DTYPE = np.float

#### Import needed modules #####
#import pysam

def makeFragmentMat(bamfile, chrom, start, end, lower, upper, atac = 1):
    nrow = upper - lower
    ncol = end - start
    mat = np.zeros( (nrow, ncol), dtype = DTYPE)
    bamHandle = AlignmentFile(bamfile)
   
    
    for read in bamHandle.fetch(chrom, max(0, start - upper), end + upper):
        if read.is_proper_pair and not read.is_reverse:
            if atac:
                #get left position
                l_pos = read.pos + 4
                #get insert size
                #correct by 8 base pairs to be inserion to insertion
                ilen = abs(read.template_length) - 8
            else:
                l_pos = read.pos
                ilen = abs(read.template_length)
            row = ilen - lower
            col = (ilen-1)//2 + l_pos - start
            if col >= 0 and col < ncol and row < nrow and row >= 0:
                mat[row, col] += 1
    bamHandle.close()
    return mat

def getInsertions(bamfile, chrom, start, end, lower, upper, atac = 1):
    npos = end - start
    mat = np.zeros(npos, dtype = DTYPE)
    bamHandle = AlignmentFile(bamfile)
    for read in bamHandle.fetch(chrom, max(0, start - upper), end + upper):
        if read.is_proper_pair and not read.is_reverse:
            if atac:
                #get left position
                l_pos = read.pos + 4
                #get insert size
                #correct by 8 base pairs to be inserion to insertion
                ilen = abs(read.template_length) - 8
            else:
                l_pos = read.pos
                ilen = abs(read.template_length)
            r_pos = l_pos + ilen - 1
            if ilen >= lower and ilen < upper:
                if l_pos >= start and l_pos < end:
                    mat[l_pos - start] += 1
                if r_pos >= start and r_pos < end:
                    mat[r_pos - start] += 1
    bamHandle.close()
    return mat

def getStrandedInsertions(bamfile, chrom, start, end, lower, upper, atac = 1):
    npos = end - start
    matplus = np.zeros(npos, dtype = DTYPE)
    matminus = np.zeros(npos, dtype = DTYPE)
    bamHandle = AlignmentFile(bamfile)
    for read in bamHandle.fetch(chrom, max(0, start - upper), end + upper):
        if read.is_proper_pair and not read.is_reverse:
            if atac:
                #get left position
                l_pos = read.pos + 4
                #get insert size
                #correct by 8 base pairs to be inserion to insertion
                ilen = abs(read.template_length) - 8
            else:
                l_pos = read.pos
                ilen = abs(read.template_length)
            r_pos = l_pos + ilen - 1
            if ilen >= lower and ilen < upper:
                if l_pos >= start and l_pos < end:
                    matplus[l_pos - start] += 1
                if r_pos >= start and r_pos < end:
                    matminus[r_pos - start] += 1
    bamHandle.close()
    return (matplus, matminus)

def getAllFragmentSizes(bamfile, lower, upper, atac = 1):
    sizes = np.zeros(upper - lower, dtype= np.float)
    # loop over samfile
    bamHandle = AlignmentFile(bamfile)
    for read in bamHandle:
          if read.is_proper_pair and not read.is_reverse:
            if atac:

                #get insert size
                #correct by 8 base pairs to be inserion to insertion
                ilen = abs(read.template_length) - 8
            else:
                ilen = abs(read.template_length)
            if ilen < upper and ilen >= lower:
                sizes[ilen - lower]+=1
    bamHandle.close()
    return sizes

def getFragmentSizesFromChunkList(chunks, bamfile, lower, upper, atac = 1):
    sizes = np.zeros(upper - lower, dtype= np.float)
    # loop over samfile
    bamHandle = AlignmentFile(bamfile)
    for chunk in chunks:
        for read in bamHandle.fetch(chunk.chrom, max(0, chunk.start - upper), chunk.end + upper):
            if read.is_proper_pair and not read.is_reverse:
                if atac:
                    #get left position
                    l_pos = read.pos + 4
                    #get insert size
                    #correct by 8 base pairs to be inserion to insertion
                    ilen = abs(read.template_length) - 8
                else:
                    l_pos = read.pos
                    ilen = abs(read.template_length)
                center = l_pos + (ilen-1) // 2
                if ilen < upper and ilen >= lower and center >= chunk.start and center < chunk.end:
                    sizes[ilen - lower] += 1
    bamHandle.close()
    return sizes

