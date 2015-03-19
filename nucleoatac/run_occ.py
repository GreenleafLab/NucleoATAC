"""
Script to make nucleosome occupancy track!

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary python modules
import matplotlib as mpl
mpl.use('PS')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import traceback
import itertools
import pysam
from pyatac.utils import shell_command,read_chrom_sizes_from_bam
from pyatac.chunk import ChunkList
from nucleoatac.Occupancy import FragmentMixDistribution, OccupancyParameters, OccChunk
from pyatac.fragmentsizes import FragmentSizes

def _occHelper(arg):
    """function to get occupancy for a set of bed regions

    """
    (chunk, params) = arg
    try:
        occ = OccChunk(chunk)
        occ.process(params)
        if params.write_peaks:
            out = (occ.getNucDist(params.flank, params.occ_calc_params.upper),
                occ.occ, [occ.peaks[i] for i in sorted(occ.peaks.keys())])
        else:
            out = (occ.getNucDist(params.flank, params.occ_calc_params.upper),
                occ.occ)
        occ.removeData()
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
        traceback.print_exc()
        print()
        raise e
    return out

def _writeOcc(track_queue, out):
    out_handle = open(out + '.occ.bedgraph','a')
    try:
        for track in iter(track_queue.get, 'STOP'):
            track.write_track(out_handle)
            track_queue.task_done()
    except Exception, e:
        print('Caught exception when writing occupancy track\n')
        traceback.print_exc()
        print()
        raise e
    out_handle.close()
    return True


def _writePeaks(pos_queue, out):
    out_handle = open(out + '.occ_peaks.bed','a')
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



def run_occ(args, bases = 5000000):
    """run occupancy calling

    """
    chrs = read_chrom_sizes_from_bam(args.bam)
    chunks = ChunkList.read(args.bed, chromDict = chrs, min_offset = 75)
    fragment_dist = FragmentMixDistribution(0, upper = args.upper)
    if args.sizes is not None:
        tmp = FragmentSizes.open(args.sizes)
        fragment_dist.fragmentsizes = FragmentSizes(0, args.upper, vals = tmp.get(0,args.upper))
    else:
        fragment_dist.getFragmentSizes(args.bam, chunks)
    fragment_dist.modelNFR()
    fragment_dist.plotFits(args.out + '.occ_fit.eps')
    fragment_dist.fragmentsizes.save(args.out + '.fragmentsizes.txt')
    params = OccupancyParameters(chunks, fragment_dist, args.upper, sep = args.nuc_sep, min_occ = args.min_occ, min_reads = args.reads,
            flank = args.flank,write_peaks = args.write_peaks, out = args.out, bam = args.bam)
    sets = chunks.split(bases = bases)
    pool1 = mp.Pool(processes = max(1,args.cores-1))
    out_handle = open(args.out + '.occ.bedgraph','w')
    out_handle.close()
    write_queue = mp.JoinableQueue()
    write_process = mp.Process(target = _writeOcc, args=(write_queue, args.out))
    write_process.start()
    if args.write_peaks:
        peaks_handle = open(args.out + '.occ_peaks.bed','w')
        peaks_handle.close()
        peaks_queue = mp.JoinableQueue()
        peaks_process = mp.Process(target = _writePeaks, args=(peaks_queue, args.out))
        peaks_process.start()
    nuc_dist = np.zeros(args.upper)
    for j in sets:
        tmp = pool1.map(_occHelper, zip(j,itertools.repeat(params)))
        for result in tmp:
            nuc_dist += result[0]
            write_queue.put(result[1])
            if args.write_peaks:
                peaks_queue.put(result[2])
    pool1.close()
    pool1.join()
    write_queue.put('STOP')
    write_process.join()
    if args.write_peaks:
        peaks_queue.put('STOP')
        peaks_process.join()
    pysam.tabix_compress(args.out + '.occ.bedgraph', args.out + '.occ.bedgraph.gz',force = True)
    shell_command('rm ' + args.out + '.occ.bedgraph')
    pysam.tabix_index(args.out + '.occ.bedgraph.gz', preset = "bed", force = True)

    dist_out = FragmentSizes(0, args.upper, vals = nuc_dist)
    dist_out.save(args.out + '.nuc_dist.txt')

    print "Making figure"
    #make figure
    fig = plt.figure()
    plt.plot(range(0,args.upper),dist_out.get(0,args.upper),label = "Nucleosome Distribution")
    plt.xlabel("Fragment Size")
    plt.ylabel("Frequency")
    fig.savefig(args.out+'.nuc_dist.eps')
    plt.close(fig)











