"""
Script to make insertion track

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
from pyatac.tracks import InsertionTrack
from pyatac.chunk import ChunkList
from pyatac.utils import shell_command, read_chrom_sizes_from_bam

def _insHelperSmooth(arg):
    """Computes smoothed insertion track for a particular set of bed regions"""
    (chunk, args) = arg
    try:
        offset = args.smooth / 2
        ins = InsertionTrack(chunk.chrom, chunk.start - offset, chunk.end + offset)
        ins.calculateInsertions( args.bam, lower = args.lower, upper = args.upper, atac = args.atac)
        ins.smooth_track(args.smooth, window = "gaussian", mode= 'valid')
    except Exception as e:
        print(('Caught exception when processing:\n'+  chunk.asBed()+"\n"))
        traceback.print_exc()
        print()
        raise e
    return ins


def _insHelper(arg):
    (chunk, args) = arg
    try:
        ins = InsertionTrack(chunk.chrom, chunk.start, chunk.end)
        ins.calculateInsertions( args.bam, lower = args.lower, upper = args.upper, atac = args.atac)
    except Exception as e:
        print(('Caught exception when processing:\n'+  chunk.asBed()+"\n"))
        traceback.print_exc()
        print()
        raise e
    return ins

def _writeIns(track_queue, out):
    out_handle = open(out + '.ins.bedgraph','a')
    try:
        for track in iter(track_queue.get, 'STOP'):
            track.write_track(out_handle)
            track_queue.task_done()
    except Exception as e:
        print('Caught exception when writing insertion track\n')
        traceback.print_exc()
        print()
        raise e
    out_handle.close()
    return True

def get_ins(args, bases = 50000, splitsize = 1000):
    """function to get insertions

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
    out_handle = open(args.out + '.ins.bedgraph','w')
    out_handle.close()
    write_queue = mp.JoinableQueue(maxsize = maxQueueSize)
    write_process = mp.Process(target = _writeIns, args=(write_queue, args.out))
    write_process.start()
    for j in sets:
        if args.smooth:
            tmp = pool1.map(_insHelperSmooth, list(zip(j,itertools.repeat(args))))
        else:
            tmp = pool1.map(_insHelper, list(zip(j,itertools.repeat(args))))
        for track in tmp:
            write_queue.put(track)
    pool1.close()
    pool1.join()
    write_queue.put('STOP')
    write_process.join()
    pysam.tabix_compress(args.out + '.ins.bedgraph', args.out + '.ins.bedgraph.gz', force = True)
    shell_command('rm ' + args.out + '.ins.bedgraph')
    pysam.tabix_index(args.out + '.ins.bedgraph.gz', preset = "bed", force = True)


