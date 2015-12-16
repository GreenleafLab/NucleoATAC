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
from pyatac.utils import shell_command, read_chrom_sizes_from_fasta, read_chrom_sizes_from_bam
from pyatac.chunk import ChunkList
from nucleoatac.NFRCalling import NFRParameters, NFRChunk
from pyatac.bias import PWM

def _nfrHelper(arg):
    """function to get occupancy for a set of bed regions

    """
    (chunk, params) = arg
    try:
        nfr = NFRChunk(chunk)
        nfr.process(params)
        if params.ins_track is None:
            out = (nfr.nfrs, nfr.ins)
        else:
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

def _writeIns(track_queue, out):
    out_handle = open(out + '.ins.bedgraph','a')
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



def run_nfr(args):
    """run nfr calling

    """
    if args.bam is None and args.ins_track is None:
        raise Exception("Must supply either bam file or insertion track")
    if not args.out:
        args.out = '.'.join(os.path.basename(args.calls).split('.')[0:-3])
    if args.fasta is not None:
        chrs_fasta = read_chrom_sizes_from_fasta(args.fasta)
        pwm = PWM.open(args.pwm)
        chunks = ChunkList.read(args.bed, chromDict = chrs_fasta, min_offset = max(pwm.up, pwm.down))
    else:
        chunks = ChunkList.read(args.bed)
    if args.bam is not None:
        chrs_bam = read_chrom_sizes_from_bam(args.bam)
        chunks.checkChroms(chrs_bam, chrom_source = "BAM file") 
    chunks.merge()
    maxQueueSize = args.cores * 10 
    params = NFRParameters(args.occ_track, args.calls, args.ins_track, args.bam, max_occ = args.max_occ, max_occ_upper = args.max_occ_upper,
                            fasta = args.fasta, pwm = args.pwm)
    sets = chunks.split(items = args.cores * 5)
    pool1 = mp.Pool(processes = max(1,args.cores-1))
    nfr_handle = open(args.out + '.nfrpos.bed','w')
    nfr_handle.close()
    nfr_queue = mp.JoinableQueue()
    nfr_process = mp.Process(target = _writeNFR, args=(nfr_queue, args.out))
    nfr_process.start()
    if params.ins_track is None:
        ins_handle = open(args.out + '.ins.bedgraph','w')
        ins_handle.close()
        ins_queue = mp.JoinableQueue()
        ins_process = mp.Process(target = _writeIns, args=(ins_queue, args.out))
        ins_process.start()
    for j in sets:
        tmp = pool1.map(_nfrHelper, zip(j,itertools.repeat(params)))
        for result in tmp:
            if params.ins_track is None:
                nfr_queue.put(result[0])
                ins_queue.put(result[1])
            else:
                nfr_queue.put(result)
    pool1.close()
    pool1.join()
    nfr_queue.put('STOP')
    nfr_process.join()
    if params.ins_track is None:
        ins_queue.put('STOP')
        ins_process.join()
    pysam.tabix_compress(args.out + '.nfrpos.bed', args.out + '.nfrpos.bed.gz',force = True)
    shell_command('rm ' + args.out + '.nfrpos.bed')
    pysam.tabix_index(args.out + '.nfrpos.bed.gz', preset = "bed", force = True)
    if params.ins_track is None:
        pysam.tabix_compress(args.out + '.ins.bedgraph', args.out + '.ins.bedgraph.gz', force = True)
        shell_command('rm ' + args.out + '.ins.bedgraph')
        pysam.tabix_index(args.out + '.ins.bedgraph.gz', preset = "bed", force = True)








