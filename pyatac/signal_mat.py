#!/usr/bin/env python
"""
Get Matrix of Sites

@author: Alicia Schep, Greenleaf lab at Stanford University
"""

##### IMPORT MODULES #####
# import necessary for python
import numpy as np
import matplotlib.pyplot as plt
from pyatac.bedgraph import BedGraphFile
from multiprocessing import Pool
from itertools import repeat
import traceback

def _signalHelperAll(arg):
    """Get aggregate signal around a set of sites

    """
    (chunks, params) = arg
    bg = BedGraphFile(params.bg)
    try:
        chunk.center()
        chunk.slop(chromDict = params.chrs, up = params.up, down = params.down)
        sig = bg.read(chunk.chrom, chunk.start, chunk.end)
        if len(sig) != (params.up + params.down + 1):
            if chunk.start == 0 :
                mat = np.hstack(np.zeros(params.up + params.down + 1 - len(sig])),
                                    sig)
            else:
                mat = np.hstack(sig,np.zeros(params.up + params.down + 1 - len(sig)))
        if chunk.strand == "-":
            sig = sig[::-1]
        if params.exp:
            sig = np.exp(sig)
        if sum(np.isnan(track))==0:
            if params.trunc:
                sig[sig<0]=0
            if params.scale:
                sig = sig / (np.sum(abs(sig))+ (np.sum(abs(sig))==0))
    except Exception as e:
        print('Caught exception when processing:\n' + chunk.asBed() + "\n")
        traceback.print_exc()
        print()
        raise e
    return sig

def _signalParams:
    def __init__(self,bg, up ,down, na):
        self.bg = bg
        self.up = up
        self.down = down
        self.na = na

def signal_mat(bed, up, down, strand_col = None, na = 0, cores = 2):
    """function to get signal from a track around some sites

    """
    chunks = read_chunks(bed, strand_col = strand_col)
    pool = Pool(processes = cores)
    result = pool.map(_signalHelperAll, zip(sets,itertools.repeat(params)))
    pool.close()
    pool.join()




