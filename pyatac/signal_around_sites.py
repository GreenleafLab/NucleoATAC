"""
Function to make plot of aggregate signal around sites.

@author: Alicia Schep, Greenleaf lab at Stanford University
"""

##### IMPORT MODULES #####
# import needed python modules
import os
import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import itertools
import traceback

#Import needed functions/classes from pyatac
from pyatac.bedgraph import BedGraphFile
from pyatac.chunk import ChunkList
from pyatac.utils import read_chrom_sizes


def _signalHelper(arg):
    """Get aggregate signal around a set of sites

    """
    (chunks, params) = arg
    bg = BedGraphFile(params.bg)
    agg = np.zeros(params.up + params.down + 1)
    if params.all:
        counter = 0
        mat = np.zeros((len(chunks),params.up + params.down + 1))
    else:
        agg = np.zeros(params.up + params.down + 1)
    for chunk in chunks:
        try:
            chunk.center()
            if params.up!=0 and params.down!=0:
                chunk.slop(chromDict = params.chrs, up = params.up, down = params.down)
            sig = bg.read(chunk.chrom, chunk.start, chunk.end)
            if params.up!=0 and params.down!=0 and len(sig) != (params.up + params.down + 1):
                if chunk.start == 0:
                    sig = np.hstack((np.zeros(params.up + params.down + 1 - len(sig)),
                                    sig))
                else:
                    sig = np.hstack((sig,np.zeros(params.up + params.down + 1 - len(sig))))
            if chunk.strand == "-":
                sig = sig[::-1]
            if params.exp:
                sig = np.exp(sig)
            if params.positive:
                sig[sig<0]=0
            if params.scale:
                tmp = sig
                tmp[np.isnan(tmp)]=0
                sig = sig / (np.sum(abs(sig))+ (np.sum(abs(sig))==0))
            if params.all:
                mat[counter] = sig
                counter += 1
            else:
                sig[np.isnan(sig)]=0
                agg += sig
        except Exception as e:
            print('Caught exception when processing:\n' + chunk.asBed() + "\n")
            traceback.print_exc()
            print()
            bg.close()
            raise e
    bg.close()
    if params.all:
        return mat
    else:
        return agg


class _signalParams:
    def __init__(self, bedgraph, sizes, up, down, exp, scale, positive, all):
        self.bg = bedgraph
        self.up = up
        self.down = down
        self.chrs = read_chrom_sizes(sizes)
        self.scale = scale
        self.positive = positive
        self.exp = exp
        self.all = all

def get_signal(args):
    """function to get signal from a track around some sites

    """
    if not args.out:
        args.out = '.'.join(os.path.basename(args.bed).split('.')[0:-1])
    chunks = ChunkList.read(args.bed, strand_col = args.strand)
    params = _signalParams(args.bg, args.sizes, args.up, args.down,args.exp,
                             args.scale, args.positive, args.all)
    sets = chunks.split(items = min(args.cores*20,len(chunks)))
    pool = Pool(processes = args.cores)
    tmp = pool.map(_signalHelper, zip(sets,itertools.repeat(params)))
    pool.close()
    pool.join()
    if args.all:
        mat = np.vstack(tmp)
        np.savetxt(args.out + ".tracks.txt.gz", mat, delimiter = ",", fmt="%1.5g")
        mat[np.isnan(mat)]=0
        result = np.sum(mat, axis = 0)
    else:
        result = sum(tmp)
    if not args.no_agg:
        if args.norm:
            result = result / len(chunks)
        fig = plt.figure()
        plt.plot(range(-args.up,args.down+1),result)
        plt.xlabel("Position relative to Site")
        plt.ylabel("Signal Intensity")
        fig.savefig(args.out+'.agg.track.eps')
        plt.close(fig)
        np.savetxt(args.out+'.agg.track.txt',result,delimiter="\t")




