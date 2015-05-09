"""
Script to determine NFR positions!

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary python modules
import multiprocessing as mp
import numpy as np
import os
import traceback
import itertools
import pysam
from pyatac.utils import shell_command
from pyatac.chunk import ChunkList
from nucleoatac.NFRCalling import NFRParameters, NFRChunk

def _nfrHelper(arg):
    """function to get occupancy for a set of bed regions

    """
    (chunk, params) = arg
    try:
        nfr = NFRChunk(chunk)
        nfr.process(params)
        out = nfr.nfrs
        nfr.removeData()
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
        traceback.print_exc()
        print()
        raise e
    return out


def _writeNFR(pos_queue, out):
    out_handle = open(out + '.nfrpos.bed','a')
    try:
        for poslist in iter(pos_queue.get, 'STOP'):
            for pos in poslist:
                pos.write(out_handle)
            pos_queue.task_done()
    except Exception, e:
        print('Caught exception when writing occupancy track\n')
        traceback.print_exc()
        print()
        raise e
    out_handle.close()
    return True



def run_nfr(args):
    """run nfr calling

    """
    if not args.out:
        args.out = '.'.join(os.path.basename(args.calls).split('.')[0:-3])
    chunks = ChunkList.read(args.bed)
    chunks.merge()
    maxQueueSize = args.cores * 10 
    params = NFRParameters(args.occ_track, args.calls, max_occ = args.max_occ, max_occ_upper = args.max_occ_upper, 
                            max_nfr_gap = args.max_nfr_gap, min_nfr_len = args.min_nfr_len)
    sets = chunks.split(items = args.cores * 5)
    pool1 = mp.Pool(processes = max(1,args.cores-1))
    nfr_handle = open(args.out + '.nfrpos.bed','w')
    nfr_handle.close()
    nfr_queue = mp.JoinableQueue()
    nfr_process = mp.Process(target = _writeNFR, args=(nfr_queue, args.out))
    nfr_process.start()
    for j in sets:
        tmp = pool1.map(_nfrHelper, zip(j,itertools.repeat(params)))
        for result in tmp:
            nfr_queue.put(result)
    pool1.close()
    pool1.join()
    nfr_queue.put('STOP')
    nfr_process.join()
    pysam.tabix_compress(args.out + '.nfrpos.bed', args.out + '.nfrpos.bed.gz',force = True)
    shell_command('rm ' + args.out + '.nfrpos.bed')
    pysam.tabix_index(args.out + '.nfrpos.bed.gz', preset = "bed", force = True)







