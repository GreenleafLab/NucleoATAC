#!/usr/bin/env python
"""
Created on Mon May 12 14:40:53 2014

@author: alicia
"""

from __future__ import division
import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.float
ctypedef np.float_t DTYPE_t



@cython.boundscheck(False)
def calculateCov(np.ndarray[DTYPE_t, ndim=1] p, np.ndarray[DTYPE_t, ndim=1] v, int r):
    if p.shape[0]!= v.shape[0]:
        raise ValueError("p and v must be same shape")
    cdef DTYPE_t value
    cdef unsigned int i, j
    for i in range(p.shape[0]):
        for j in range(i,p.shape[0]):
            if i==j:
                value += p[i]* (1-p[i]) * v[i]**2
            else:
                value += p[i] * p[j] * -2 * v[i] * v[j]
    return value * r
