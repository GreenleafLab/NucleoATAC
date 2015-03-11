#!/usr/bin/env python
"""
Script to call nucleosome positions-- track making, nucleosome calling, and nfr calling!

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary python modules
import matplotlib as mpl
mpl.use('PS')
import PyATAC as PA
import VMat as V
import NucleosomeCalling as Nuc
from multiprocessing import Pool
import Occupancy as Occ
import traceback
from itertools import repeat

def nuc_helper(arg):
    """function to call nucs for particular set of bed regions

    """
    (chunkset, params) = arg
    try:
        chunkset.process(params)
    except Exception as e:
        print('Caught exception when processing: '+  chunkset.name)
        traceback.print_exc()
        print()
        raise e
 


def run_nuc(args):
    """run nucleosome calling-- make signal tracks, call nucleosome positions

    """
    #Make tmp directory for tmp files
    try:
        PA.shell_command('mkdir tmp_nucleoatac')
    except:
        print "Deleting old tmp_nuc folder to make new one"
        PA.shell_command('rm -r tmp_nucleoatac')
        PA.shell_command('mkdir tmp_nucleoatac')
    
    print "Reading in bedfile and setting up parameters"
    #Set up parameters
    ##Get chromosome info from bam
    chromosomes = PA.Chromosomes()
    chromosomes.get_from_bam(args.bam)
    chromosomes.write_chrs('tmp_nucleoatac/genome.txt')
    
    vmat = V.VMat.open(args.vmat)
    
    bed_list = Nuc.read_bed(args.bed, chromosomes.chrs,
                        max_offset = vmat.mat.shape[1] + vmat.i_upper/2)
    sets = Nuc.makeChunkSetList(bed_list, args.cores)
    
    insert_dist = Occ.InsertDistribution(0,i_upper = vmat.i_upper)
    insert_dist.getInsertSizes(args.bam, bed_list)

    if args.occ_track is None:
        insert_dist.modelNFR(vmat.i_lower)
        occ_params = Occ.OccupancyCalcParams(0, vmat.i_upper, insert_dist)
    else:
        occ_params = None
    
    params = Nuc.NucParameters(vmat, insert_dist.insertsizes, sd = args.sd, nonredundant_sep = args.nuc_sep,
                 redundant_sep = args.redundant_sep, min_z = args.min_z, min_lr = args.min_lr, downsample = args.downsample,
                 seed = args.seed, min_nfr_len = 1, max_nfr_len = 1000, max_nfr_occ = args.max_nfr_occ,
                out = args.out, gdna = args.gdna, bias= args.bias, occ_track = args.occ_track, bam = args.bam,
                write_all = args.write_all, occ_params = occ_params)
    
    print "Processing Regions"
    
    ##parallel processed computation
    pool = Pool(processes=args.cores)
    pool.map(nuc_helper,zip(sets,repeat(params)))
    
    print "Merging tmp files"
    
    if args.write_all:
        track_outputs = ['nucsig','smoothsig','ins','nuc_cov','nfr_cov','background','rawsig','fitted']
    else:
        track_outputs = ['nucsig','smoothsig','ins']
    if args.occ_track is None:
        track_outputs.append('occ')
    
    for track in track_outputs:
        PA.combine_tmp_bed_to_bw(args.out, track)
    
    for bed in ['nucpos','nucpos.redundant','nfrpos']:
        PA.combine_tmp_bed(args.out, bed)
    
    PA.shell_command('rm -r tmp_nucleoatac')









