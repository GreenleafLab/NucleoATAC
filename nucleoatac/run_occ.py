#!/usr/bin/env python
"""
Script to make nucleosome occupancy track!

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary python modules
import matplotlib as mpl
mpl.use('PS')
import matplotlib.pyplot as plt
import PyATAC as PA
import Occupancy as Occ
from multiprocessing import Pool
import numpy as np
from itertools import repeat
import traceback

def occ_helper(arg):
    """function to get occupancy for a set of bed regions

    """
    (chunkset, params) = arg
    try:
        nuc_dist = chunkset.process(params)
    except Exception as e:
        print('Caught exception when processing: '+  chunkset.name)
        traceback.print_exc()
        print()
        raise e
    return nuc_dist

def run_occ(args):
    """run occupancy calling 

    """
    #Make tmp directory for tmp files
    try:
        PA.shell_command('mkdir tmp_nucleoatac')
    except:
        print "Deleting old tmp_nuc folder to make new one"
        PA.shell_command('rm -r tmp_nucleoatac')
        PA.shell_command('mkdir tmp_nucleoatac')
    
    
    ##Get chromosome info from bam
    chromosomes = PA.Chromosomes()
    chromosomes.get_from_bam(args.bam)
    chromosomes.write_chrs('tmp_nucleoatac/genome.txt')
    
    #Read in bed list
    bed_list = Occ.read_bed(args.bed, chromosomes.chrs, 75)
    
    print "Getting insert distribution"
    insert_dist = Occ.InsertDistribution(0,i_upper = args.upper)
    insert_dist.getInsertSizes(args.bam, bed_list)
    
    print "Modelling NFR distribution"
    insert_dist.modelNFR(args.lower)
    insert_dist.plotFits(args.out+".fits.eps")
    
    print "Calculating occupancies and getting nucleosomal distribution"
    
    params = Occ.OccupancyParams(Occ.OccupancyCalcParams(0,args.upper, insert_dist),
                            sep = args.nuc_sep, min_occ = args.min_occ, min_reads = args.reads, flank = args.flank, 
                            write_peaks = args.write_peaks, out = args.out, bam = args.bam)
    
    #divide intervals into sets based on number of cores
    sets = Occ.makeChunkSetList(bed_list, args.cores)

    ##parallel processed computation
    pool = Pool(processes=args.cores)
    submats=pool.map(occ_helper,zip(sets,repeat(params)))
    
    mat = np.zeros(submats[0].shape)
    for i in range(len(submats)):
        mat=mat+submats[i]
    mat = mat / sum(mat)
    nuc_dist = PA.InsertSizes(mat,0,args.upper)
    nuc_dist.save(args.out + '.nuc_dist.txt')
    
    print "Making figure"
    #make figure
    fig = plt.figure()
    plt.plot(range(0,args.upper),nuc_dist.get(0,args.upper),label = "Nucleosome Distribution")
    plt.xlabel("Fragment Size")
    plt.ylabel("Frequency")
    fig.savefig(args.out+'.nuc_dist.eps')
    plt.close(fig)
    
    print "Merging tmp files"
    
    PA.combine_tmp_bed_to_bw(args.out, 'occ')
    
    if args.write_peaks:
        PA.combine_tmp_bed(args.out, 'occ_peaks')

    PA.shell_command('rm -r tmp_nucleoatac')










