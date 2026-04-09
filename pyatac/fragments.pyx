# cython: language_level=3

import pyatac.seq as seq
import numpy as np
cimport numpy as np
cimport cython
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


@cython.boundscheck(False)
def makeFragmentMat(str bamfile, str chrom, int start, int end,
                    int lower, int upper, int atac=1):

    cdef int nrow = upper - lower
    cdef int ncol = end - start
    cdef np.ndarray[DTYPE_t, ndim=2] mat = np.zeros((nrow, ncol), dtype=DTYPE)

    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int l_pos, ilen, row, col

    for read in bamHandle.fetch(chrom, max(0, start - upper), end + upper):
        if read.is_proper_pair and not read.is_reverse:
            if atac:
                l_pos = read.reference_start + 4
                ilen = abs(read.template_length) - 8
            else:
                l_pos = read.reference_start
                ilen = abs(read.template_length)

            row = ilen - lower
            col = (ilen - 1) // 2 + l_pos - start

            if 0 <= col < ncol and 0 <= row < nrow:
                mat[row, col] += 1

    bamHandle.close()
    return mat


@cython.boundscheck(False)
def getInsertions(str bamfile, str chrom, int start, int end,
                  int lower, int upper, int atac=1):

    cdef int npos = end - start
    cdef np.ndarray[DTYPE_t, ndim=1] mat = np.zeros(npos, dtype=DTYPE)

    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int l_pos, ilen, r_pos

    for read in bamHandle.fetch(chrom, max(0, start - upper), end + upper):
        if read.is_proper_pair and not read.is_reverse:
            if atac:
                l_pos = read.reference_start + 4
                ilen = abs(read.template_length) - 8
            else:
                l_pos = read.reference_start
                ilen = abs(read.template_length)

            r_pos = l_pos + ilen - 1

            if lower <= ilen < upper:
                if start <= l_pos < end:
                    mat[l_pos - start] += 1
                if start <= r_pos < end:
                    mat[r_pos - start] += 1

    bamHandle.close()
    return mat


@cython.boundscheck(False)
def getStrandedInsertions(str bamfile, str chrom, int start, int end,
                          int lower, int upper, int atac=1):

    cdef int npos = end - start
    cdef np.ndarray[DTYPE_t, ndim=1] matplus = np.zeros(npos, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] matminus = np.zeros(npos, dtype=DTYPE)

    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int l_pos, ilen, r_pos

    for read in bamHandle.fetch(chrom, max(0, start - upper), end + upper):
        if read.is_proper_pair and not read.is_reverse:
            if atac:
                l_pos = read.reference_start + 4
                ilen = abs(read.template_length) - 8
            else:
                l_pos = read.reference_start
                ilen = abs(read.template_length)

            r_pos = l_pos + ilen - 1

            if lower <= ilen < upper:
                if start <= l_pos < end:
                    matplus[l_pos - start] += 1
                if start <= r_pos < end:
                    matminus[r_pos - start] += 1

    bamHandle.close()
    return (matplus, matminus)


@cython.boundscheck(False)
def getAllFragmentSizes(str bamfile, int lower, int upper, int atac=1):

    cdef np.ndarray[DTYPE_t, ndim=1] sizes = np.zeros(upper - lower, dtype=DTYPE)

    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int ilen

    for read in bamHandle:
        if read.is_proper_pair and not read.is_reverse:
            if atac:
                ilen = abs(read.template_length) - 8
            else:
                ilen = abs(read.template_length)

            if lower <= ilen < upper:
                sizes[ilen - lower] += 1

    bamHandle.close()
    return sizes


@cython.boundscheck(False)
def getFragmentSizesFromChunkList(chunks, str bamfile,
                                  int lower, int upper, int atac=1):

    cdef np.ndarray[DTYPE_t, ndim=1] sizes = np.zeros(upper - lower, dtype=DTYPE)

    cdef AlignmentFile bamHandle = AlignmentFile(bamfile)
    cdef AlignedSegment read
    cdef int ilen, l_pos, center

    for chunk in chunks:
        for read in bamHandle.fetch(chunk.chrom,
                                    max(0, chunk.start - upper),
                                    chunk.end + upper):

            if read.is_proper_pair and not read.is_reverse:
                if atac:
                    l_pos = read.reference_start + 4
                    ilen = abs(read.template_length) - 8
                else:
                    l_pos = read.reference_start
                    ilen = abs(read.template_length)

                center = l_pos + (ilen - 1) // 2

                if (lower <= ilen < upper and
                        chunk.start <= center < chunk.end):
                    sizes[ilen - lower] += 1

    bamHandle.close()
    return sizes
