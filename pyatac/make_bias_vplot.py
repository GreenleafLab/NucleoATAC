#!/usr/bin/env python
"""
Script to make V-plot

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary python modules
import numpy as np
import os
#import matplotlib as mpl
#mpl.use('PS')
from multiprocessing import Pool
import traceback
import itertools
#import pyatac functions
from pyatac.utils import read_chrom_sizes_from_fasta
from pyatac.bias import InsertionBiasTrack, PWM
from pyatac.chunkmat2d import BiasMat2D
from pyatac.chunk import ChunkList
from pyatac.VMat import VMat
from pyatac.fragmentsizes import FragmentSizes

def _biasVplotHelper(arg):
    """function to make vplot for a particular set of bed regions

    """
    (chunks, params) = arg
    mat = np.zeros((params.upper-params.lower,2*params.flank+1))
    for chunk in chunks:
        try:
            chunk.center()
            biastrack = InsertionBiasTrack(chunk.chrom, chunk.start - params.flank - 1 - (params.upper/2),
                                                chunk.end + params.flank + params.upper/2+1)
            if params.bg is not None:
                biastrack.read_track(params.bg, empty = 0)
            else:
                biastrack.computeBias(params.fasta, params.chrs, params.pwm)
            biasmat = BiasMat2D(chunk.chrom, chunk.start - params.flank - 1, chunk.end + params.flank,
                                                params.lower, params.upper)
            biasmat.makeBiasMat(biastrack)
            biasmat.normByInsertDist(params.fragmentsizes)
            add = biasmat.get(start = chunk.start - params.flank, end = chunk.end + params.flank,
                                flip = (chunk.strand == "-"))
            if params.scale:
                mat += add/np.sum(add)
            else:
                mat += add
        except Exception as e:
            print('Caught exception when processing:\n' + chunk.asBed() + "\n")
            traceback.print_exc()
            print()
            raise e
    return mat

class _BiasVplotParams:
    def __init__(self, flank, lower, upper, bg, fasta, pwm, sizes, scale):
        self.flank = flank
        self.lower = lower
        self.upper = upper
        self.scale = scale
        self.bg = bg
        self.fasta = fasta
        if self.bg is None:
            self.pwm = PWM.open(pwm)
            self.chrs = read_chrom_sizes_from_fasta(fasta)
        self.fragmentsizes = FragmentSizes.open(sizes)


def make_bias_vplot(args):
    """function to make vplot

    """
    if not args.out:
        args.out = '.'.join(os.path.basename(args.bed).split('.')[0:-1])
    chunks = ChunkList.read(args.bed, strand_col = args.strand)
    sets = chunks.split(items = min(args.cores*20,len(chunks)))
    params = _BiasVplotParams(flank = args.flank, lower = args.lower, upper = args.upper, bg = args.bg,
                                sizes = args.sizes, scale = args.scale,
                                pwm = args.pwm, fasta = args.fasta)
    pool = Pool(processes = args.cores)
    tmp = pool.map(_biasVplotHelper, zip(sets,itertools.repeat(params)))
    pool.close()
    pool.join()
    result = sum(tmp)
    ##Turn matrix into VMat object
    vmat=VMat(result,args.lower,args.upper)
    vmat.plot(filename=args.out+".Bias.Vplot.eps")
    if args.plot_extra:
        ##get insertion profile represented by vplot
        vmat.converto1d()
        vmat.plot_1d(filename=args.out+'.InsertionProfile.eps')
        #get insert size dstribution represented by vplot
        vmat.plot_insertsize(filename= args.out + ".InsertSizes.eps")
    ##save
    vmat.save(args.out+".Bias.VMat")


