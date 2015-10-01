"""
Script to compute Tn5 bias

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary for python
import os
import multiprocessing as mp
import itertools
import traceback
import pysam
import numpy as np
from pyatac.chunk import ChunkList
from pyatac.bias import InsertionBiasTrack, PWM
from pyatac.utils import read_chrom_sizes_from_fasta, shell_command

def _biasHelper(arg):
    """Helper function to multiprocess computation of bias tracks"""
    (chunk, params) = arg
    try:
        bias = InsertionBiasTrack(chunk.chrom, chunk.start, chunk.end)
        bias.computeBias(params.fasta, params.chrs, params.pwm)
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
        traceback.print_exc()
        print()
        raise e
    return bias

def _writeBias(track_queue, out):
    """Function to handle writing of bias output"""
    out_handle = open(out + '.Scores.bedgraph','a')
    try:
        for track in iter(track_queue.get, 'STOP'):
            track.write_track(out_handle)
            track_queue.task_done()
    except Exception, e:
        print('Caught exception when writing insertion track\n')
        traceback.print_exc()
        print()
        raise e
    out_handle.close()
    return True

class _BiasParams:
    """Class to store parameter for _biasHelper function"""
    def __init__(self, fasta, pwm):
        self.chrs = read_chrom_sizes_from_fasta(fasta)
        self.fasta = fasta
        self.pwm = PWM.open(pwm)

def make_bias_track(args, bases = 500000, splitsize = 1000):
    """function to compute bias track

    """
    if args.out is None:
        if args.bed is not None:
            args.out = '.'.join(os.path.basename(args.bed).split('.')[0:-1])
        else:
            args.out = '.'.join(os.path.basename(args.fasta).split('.')[0:-1])
    params = _BiasParams(args.fasta, args.pwm)
    if args.bed is None:
        chunks = ChunkList.convertChromSizes(params.chrs, splitsize = splitsize)
        sets = chunks.split(items = bases/splitsize)
    else:
        chunks = ChunkList.read(args.bed)
        chunks.checkChroms(params.chrs.keys())
        chunks.merge()
        sets = chunks.split(bases = bases)
    maxQueueSize = max(2,int(2 * bases / np.mean([chunk.length() for chunk in chunks])))
    pool = mp.Pool(processes = max(1,args.cores-1))
    out_handle = open(args.out + '.Scores.bedgraph','w')
    out_handle.close()
    write_queue = mp.JoinableQueue(maxsize = maxQueueSize)
    write_process = mp.Process(target = _writeBias, args=(write_queue, args.out))
    write_process.start()
    for j in sets:
        tmp = pool.map(_biasHelper, zip(j,itertools.repeat(params)))
        for track in tmp:
            write_queue.put(track)
    pool.close()
    pool.join()
    write_queue.put('STOP')
    write_process.join()
    pysam.tabix_compress(args.out + '.Scores.bedgraph', args.out + '.Scores.bedgraph.gz', force = True)
    shell_command('rm ' + args.out + '.Scores.bedgraph')
    pysam.tabix_index(args.out + '.Scores.bedgraph.gz', preset = "bed", force = True)

