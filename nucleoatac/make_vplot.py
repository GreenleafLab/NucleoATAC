#!/usr/bin/env python
"""
Script to make V-plot

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary for python
import numpy as np
import matplotlib as mpl
mpl.use('PS')
import pysam
import PyATAC as PA
import VMat as V
from multiprocessing import Pool
from itertools import repeat
import traceback


def vplot_helper(arg):
    """function to make vplot for a particular set of bed regions

    """
    (chunkset, args) = arg
    try:
        bamfile = pysam.Samfile(args.bam, "rb")
        mat = np.zeros((args.upper-args.lower,2*args.flank+1))
        for chunk in chunkset.chunks:
            if chunk.length > 1:
                mid=(chunk.start+(chunk.end-chunk.start-1)/2)
                chunk.start = mid
                chunk.end = mid + 1
                chunk.length = 1
            elif chunk.length == 0:
                chunk.end = chunk.start + 1
                chunk.length = 1
            reads = PA.ReadList(chunk.chrom, chunk.start - args.flank-1, chunk.end + args.flank)
            reads.extract_reads(bamfile, i_upper = args.upper, atac = args.not_atac)
            if chunk.strand=="-":
                submat = PA.ReadMat2D(chunk.chrom, chunk.start - args.flank-1, chunk.end + args.flank, args.lower, args.upper)
            else:
                submat = PA.ReadMat2D(chunk.chrom, chunk.start - args.flank, chunk.end + args.flank , args.lower, args.upper)
            submat.makeReadMat(reads)
            if np.sum(submat.mat)==0:
                continue
            if chunk.strand == "-":
                add = submat.get(start = chunk.start- args.flank, end = chunk.end + args.flank, flip = True)
            else:
                add= submat.get()
            if args.scale:
                mat += add/np.sum(add)
            else:
                mat += add
        bamfile.close()
    except Exception as e:
        print('Caught exception when processing: '+  chunkset.name)
        traceback.print_exc()
        print()
        raise e
    return mat

def make_vplot(args):
    """function to make vplot

    """
    #read in bed list
    bed_list = PA.read_bed(args.bed,weight_col=args.weight,strand_col = args.strand)
    #divide intervals into sets based on number of cores
    sets = PA.makeChunkSetList(bed_list, args.cores)

    ##parallel processed computation
    pool = Pool(processes=args.cores)
    submats=pool.map(vplot_helper,zip(sets,repeat(args)))
    # sum up matrices for each chunk into matrix for all
    mat = np.zeros(submats[0].shape)
    for i in range(len(submats)):
        mat=mat+submats[i]

    ##Turn matrix into VMat object
    vmat=V.VMat(mat,args.lower,args.upper)
    vmat.plot(filename=args.out+".Vplot.eps")
    ##get insertion profile represented by vplot
    vmat.converto1d()
    vmat.plot_1d(filename=args.out+'.InsertionProfile.eps')
    #get insert size dstribution represented by vplot
    vmat.plot_insertsize(filename= args.out + ".InsertSizes.eps")
    ##save
    vmat.save(args.out+".VMat")











