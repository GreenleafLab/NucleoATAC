#### Import needed modules #####
import pyatac.seq as seq
import numpy as np
cimport numpy as np
cimport cython
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam.libcfaidx cimport FastaFile

DTYPE = np.float
ctypedef np.float_t DTYPE_t


#### Import needed modules #####
#import pysam

@cython.boundscheck(False)
def makeFragmentMat(str bamfile, str chrom, int start, int end, int lower, int upper, int atac = 1):
    cdef int nrow = upper - lower
    cdef int ncol = end - start
    cdef np.ndarray[DTYPE_t, ndim=2] mat = np.zeros( (nrow, ncol), dtype = DTYPE)
    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int l_pos, ilen, row, col
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
            col = (ilen-1)/2 + l_pos - start
            if col >= 0 and col < ncol and row < nrow and row >= 0:
                mat[row, col] += 1
    bamHandle.close()
    return mat

@cython.boundscheck(False)
def getInsertions(str bamfile, str chrom, int start, int end, int lower, int upper, int atac = 1):
    cdef int npos = end - start
    cdef np.ndarray[DTYPE_t, ndim=1] mat = np.zeros(npos, dtype = DTYPE)
    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int l_pos, ilen, r_pos
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


@cython.boundscheck(False)
def getStrandedInsertions(str bamfile, str chrom, int start, int end, int lower, int upper, int atac = 1):
    cdef int npos = end - start
    cdef np.ndarray[DTYPE_t, ndim=1] matplus = np.zeros(npos, dtype = DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] matminus = np.zeros(npos, dtype = DTYPE)
    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int l_pos, ilen, r_pos
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



@cython.boundscheck(False)
def getAllFragmentSizes(str bamfile, int lower, int upper, int atac = 1):
    cdef np.ndarray[DTYPE_t, ndim =1] sizes = np.zeros(upper - lower, dtype= np.float)
    # loop over samfile
    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int ilen
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


@cython.boundscheck(False)
def getFragmentSizesFromChunkList(chunks, str bamfile, int lower, int upper, int atac = 1):
    cdef np.ndarray[DTYPE_t, ndim =1] sizes = np.zeros(upper - lower, dtype= np.float)
    # loop over samfile
    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int ilen, l_pos, center
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
                center = l_pos + (ilen-1)/2
                if ilen < upper and ilen >= lower and center >= chunk.start and center < chunk.end:
                    sizes[ilen - lower]+=1
    bamHandle.close()
    return sizes

