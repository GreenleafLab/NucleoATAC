"""
Script to make nucleosome occupancy track!

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary python modules
import multiprocessing as mp
import numpy as np
import traceback
import itertools
import pysam
from pyatac.utils import shell_command,read_chrom_sizes_from_bam
from pyatac.chunk import ChunkList
from nucleoatac.Occupancy import FragmentMixDistribution, OccupancyParameters, OccChunk
from pyatac.fragmentsizes import FragmentSizes
from pyatac.bias import PWM

def _diffHelper(arg):
    """function to get occupancy for a set of bed regions

    """
    (chunk, params) = arg
    try:
        occ = OccChunk(chunk)
        occ.process(params)
        out = (occ.getNucDist(),
                occ.occ, [occ.peaks[i] for i in sorted(occ.peaks.keys())])
        occ.removeData()
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
        traceback.print_exc()
        print()
        raise e
    return out


def _writeDiff(pos_queue, out):
    out_handle = open(out + '.occdiff.bed','a')
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

def run_diff(args, bases = 500000):
    """run differential occupancy calling

    """
    chrs = read_chrom_sizes_from_bam(args.bam)
    pwm = PWM.open(args.pwm)
    chunks = ChunkList.read(args.bed, chromDict = chrs, min_offset = args.flank + args.upper/2 + max(pwm.up,pwm.down))
    chunks.merge()
    maxQueueSize = max(2,int(100 * bases / np.mean([chunk.length() for chunk in chunks])))
    #get fragmentsizes
    fragment_dist1 = FragmentMixDistribution(0, upper = args.upper)
    fragment_dist1.fragmentsizes = FragmentSizes(0, args.upper, vals = FragmentSizes.open(args.sizes1).get(0,args.upper))
    fragment_dist1.modelNFR()
    fragment_dist2 = FragmentMixDistribution(0, upper = args.upper)
    fragment_dist2.fragmentsizes = FragmentSizes(0, args.upper, vals = FragmentSizes.open(args.sizes2).get(0,args.upper))
    fragment_dist2.modelNFR()
    params = OccupancyParameters(fragment_dist, args.upper, args.fasta, args.pwm, sep = args.nuc_sep, min_occ = args.min_occ,
            flank = args.flank, bam = args.bam, ci = args.confidence_interval)
    sets = chunks.split(bases = bases)
    pool1 = mp.Pool(processes = max(1,args.cores-1))
    diff_handle = open(args.out + '.occdiff.bed','w')
    diff_handle.close()
    diff_queue = mp.JoinableQueue()
    diff_process = mp.Process(target = _writeDiff, args=(diff_queue, args.out))
    diff_process.start()
    nuc_dist = np.zeros(args.upper)
    for j in sets:
        tmp = pool1.map(_occHelper, zip(j,itertools.repeat(params)))
        for result in tmp:
            diff_queue.put(result[1])
    pool1.close()
    pool1.join()
    diff_queue.put('STOP')
    diff_process.join()
    pysam.tabix_compress(args.out + '.occdiff.bed', args.out + '.occdiff.bed.gz',force = True)
    shell_command('rm ' + args.out + '.occdiff.bed')
    pysam.tabix_index(args.out + '.occdiff.bed.gz', preset = "bed", force = True)












