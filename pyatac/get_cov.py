"""
Script to make coverage track

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary for python
import os
import multiprocessing as mp
import itertools
import numpy as np
import pysam
import traceback
import pyximport; pyximport.install()
from pyatac.tracks import CoverageTrack
from pyatac.chunk import ChunkList
from pyatac.utils import shell_command, read_chrom_sizes_from_bam
from pyatac.chunkmat2d import FragmentMat2D

def _covHelper(arg):
    """Computes coverage track for a particular set of bed regions"""
    (chunk, args) = arg
    try:
        offset = args.window / 2
        mat = FragmentMat2D(chunk.chrom,chunk.start - offset, chunk.end + offset, args.lower, args.upper, args.atac) 
        mat.makeFragmentMat(args.bam)
        cov = CoverageTrack(chunk.chrom, chunk.start, chunk.end)
        cov.calculateCoverage(mat, lower = args.lower, upper = args.upper, window_len = args.window)
        cov.vals *= args.scale / float(args.window)
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
        traceback.print_exc()
        print()
        raise e
    return cov


def _writeCov(track_queue, out):
    out_handle = open(out + '.cov.bedgraph','a')
    try:
        for track in iter(track_queue.get, 'STOP'):
            track.write_track(out_handle)
            track_queue.task_done()
    except Exception, e:
        print('Caught exception when writing coverage track\n')
        traceback.print_exc()
        print()
        raise e
    out_handle.close()
    return True

def get_cov(args, bases = 50000, splitsize = 1000):
    """function to get coverages

    """
    if not args.out:
        if args.bed is None:
            args.out = '.'.join(os.path.basename(args.bam).split('.')[0:-1])
        else:
            args.out = '.'.join(os.path.basename(args.bed).split('.')[0:-1])
    if args.bed is None:
        chrs = read_chrom_sizes_from_bam(args.bam)
        chunks = ChunkList.convertChromSizes(chrs, splitsize = splitsize)
        sets = chunks.split(items = bases/splitsize)
    else:
        chunks = ChunkList.read(args.bed)
        chunks.merge()
        sets = chunks.split(bases = bases)
    maxQueueSize = max(2,int(2 * bases / np.mean([chunk.length() for chunk in chunks])))
    pool1 = mp.Pool(processes = max(1,args.cores-1))
    out_handle = open(args.out + '.cov.bedgraph','w')
    out_handle.close()
    write_queue = mp.JoinableQueue(maxsize = maxQueueSize)
    write_process = mp.Process(target = _writeCov, args=(write_queue, args.out))
    write_process.start()
    for j in sets:
        tmp = pool1.map(_covHelper, zip(j,itertools.repeat(args)))
        for track in tmp:
            write_queue.put(track)
    pool1.close()
    pool1.join()
    write_queue.put('STOP')
    write_process.join()
    pysam.tabix_compress(args.out + '.cov.bedgraph', args.out + '.cov.bedgraph.gz', force = True)
    shell_command('rm ' + args.out + '.cov.bedgraph')
    pysam.tabix_index(args.out + '.cov.bedgraph.gz', preset = "bed", force = True)


