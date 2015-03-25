"""
Gets aggregate nuc frequency around sites.

@author: Alicia Schep, Greenleaf lab at Stanford University
"""

##### IMPORT MODULES #####
import os
import numpy as np
import itertools
import traceback
from multiprocessing import Pool
import pyatac.seq as seq
from pyatac.chunk import ChunkList
import pyatac.utils as utils


###### Helper functions and classes ######
def _nucleotideHelper(arg):
    """Helper function for multiprocessing acquisition of sequence content around sites"""
    (chunks, params) = arg
    mat = np.zeros(params.matsize)
    n = 0.0
    try:
        for chunk in chunks:
            chunk.center()
            chunk.slop(chromDict = params.chrs, up = params.up, down = params.down + params.dinucleotide)
            sequence = seq.get_sequence(chunk, params.fasta)
            submat = seq.seq_to_mat(sequence, params.nucleotides)
            if len(sequence) == (params.up + params.down + 1 + params.dinucleotide):
                mat += submat
                n += 1
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+'\n')
        traceback.print_exc()
        print()
        raise e
    return mat,n


class _NucleotideParameters:
    """Class to store parameters related to getting nucleotides"""
    def __init__(self, up, down, fasta, dinucleotide = False):
        self.up = up
        self.down = down
        self.fasta = fasta
        self.chrs = utils.read_chrom_sizes_from_fasta(fasta)
        if dinucleotide:
            nucs="CGAT"
            dinucs = []
            for p in itertools.product(nucs, repeat=2):
                dinucs.append("".join(p))
            self.nucleotides = dinucs
        else:
            self.nucleotides = ["A","C","G","T"]
        self.matsize = (len(self.nucleotides), self.up + self.down + 1)
        self.dinucleotide = dinucleotide

##### Main function #####

def get_nucleotide(args):
    """Function to obain sequence content around sites"""
    if not args.out:
        args.out = '.'.join(os.path.basename(args.bed).split('.')[0:-1])
    chunks = ChunkList.read(args.bed, strand_col = args.strand)
    params = _NucleotideParameters(args.up, args.down, args.fasta, args.dinucleotide)
    sets = chunks.split(bases = 10000)
    pool = Pool(processes = args.cores)
    tmp = pool.map(_nucleotideHelper, zip(sets,itertools.repeat(params)))
    pool.close()
    pool.join()
    result = np.zeros(params.matsize)
    n = 0.0
    for i in tmp:
        result += i[0]
        n += i[1]
    result = result / n
    if args.norm:
        normfreqs = seq.getNucFreqs(params.fasta, params.nucleotides)
        result = result / np.reshape(np.repeat(normfreqs,result.shape[1]),result.shape)
    #save text output
    out = np.hstack((np.array(params.nucleotides)[:,np.newaxis], result.astype('|S8')))
    np.savetxt(args.out+'.nucfreq.txt', out, delimiter="\t", fmt = '%s')#,fmt = '%1.4g')



