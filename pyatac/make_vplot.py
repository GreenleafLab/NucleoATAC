"""
Script to make V-plot

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary for python
import os
import numpy as np
#import matplotlib as mpl
#mpl.use('PS')
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
from pyatac.chunk import ChunkList
from pyatac.chunkmat2d import FragmentMat2D
import VMat as V
from multiprocessing import Pool
import itertools
import traceback


def _vplotHelper(arg):
    """function to make vplot for one region

    """
    (chunks, params) = arg
    result = np.zeros(((params.upper - params.lower), 2 * params.flank + 1))
    for chunk in chunks:
        try:
            chunk.center()
            submat = FragmentMat2D(chunk.chrom, chunk.start - params.flank-1, chunk.end + params.flank, params.lower, params.upper, params.atac)
            submat.makeFragmentMat(params.bam)
            add = submat.get(start = chunk.start- params.flank, end = chunk.end + params.flank, flip = (chunk.strand =="-"))
            if params.scale:
                add = add/np.sum(add)
            result += add
        except Exception as e:
            print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
            traceback.print_exc()
            print()
            raise e
    return result


class _VplotParams:
    """Class to store parameters for use in _vplotHelper"""
    def __init__(self, flank, lower, upper, bam, atac, scale):
        self.flank = flank
        self.lower = lower
        self.upper = upper
        self.bam = bam
        self.atac = atac
        self.scale = scale





def make_vplot(args):
    """function to make vplot

    """
    if not args.out:
        args.out = '.'.join(os.path.basename(args.bed).split('.')[0:-1])
    chunks = ChunkList.read(args.bed, strand_col = args.strand)
    sets = chunks.split(items = min(args.cores*20,len(chunks)))
    params = _VplotParams(flank = args.flank, lower = args.lower, upper = args.upper, bam = args.bam,
                            atac = args.atac, scale = args.scale)
    pool = Pool(processes = args.cores)
    tmp = pool.map(_vplotHelper, zip(sets,itertools.repeat(params)))
    pool.close()
    pool.join()
    result = sum(tmp)
    ##Turn matrix into VMat object
    vmat=V.VMat(result,args.lower,args.upper)
    if not args.no_plot:
        vmat.plot(filename=args.out+".Vplot.eps")
    if args.plot_extra:
        ##get insertion profile represented by vplot
        vmat.converto1d()
        vmat.plot_1d(filename=args.out+'.InsertionProfile.eps')
        #get insert size dstribution represented by vplot
        vmat.plot_insertsize(filename= args.out + ".InsertSizes.eps")
    ##save
    vmat.save(args.out+".VMat")











