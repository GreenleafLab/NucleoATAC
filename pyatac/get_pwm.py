"""
Gets aggregate nuc frequency around sites.

@author: Alicia Schep, Greenleaf lab at Stanford University
"""

##### IMPORT MODULES #####
import os
import numpy as np
import itertools
from multiprocessing import Pool
import traceback
import pyatac.seq as seq
from pyatac.chunk import ChunkList
from pyatac.bias import PWM
import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()})
from pyatac.tracks import InsertionTrack
from pyatac.utils import read_chrom_sizes_from_fasta

def _pwmHelper(arg):
    """"Helper function for multiprocessing acquisition of sequence around insertions"""
    (chunks, params) = arg
    try:
        n = 0.0
        mat = np.zeros((len(params.nucleotides), params.up + params.down + 1))
        for chunk in chunks:
            ins = InsertionTrack(chunk.chrom, chunk.start, chunk.end)
            if params.sym:
                ins.calculateInsertions( params.bam, lower = params.lower, upper = params.upper, atac = params.atac)
                mat += ins.getInsertionSequences(params.fasta, params.nucleotides, up = params.up, down = params.down)
            else:
                ins.calculateStrandedInsertions( params.bam, lower = params.lower, upper = params.upper, atac = params.atac)
                mat += ins.getStrandedInsertionSequences(params.fasta, params.nucleotides, up = params.up, down = params.down)
            n += sum(ins.vals)
    except Exception as e:
        print('Caught exception when processing:\n'+  chunk.asBed()+"\n")
        traceback.print_exc()
        print()
        raise e
    return mat, n


class _PWMParameters:
    """Class to store parameters related to getting pwm"""
    def __init__(self, bam, up, down, fasta, lower = 0, upper = 2000, atac = True, sym = True):
        self.bam = bam
        self.up = up
        self.down = down
        self.lower = lower
        self.upper = upper
        self.fasta = fasta
        self.atac = atac
        self.nucleotides = ["A","C","G","T"]
        self.sym = sym

##### Main function #####
def get_pwm(args, bases = 50000, splitsize = 1000):
    """Functiono obtain PWM around ATAC insertion"""
    if not args.out:
        args.out = '.'.join(os.path.basename(args.bam).split('.')[0:-1])
    chrs = read_chrom_sizes_from_fasta(args.fasta)
    if args.bed is None:
        chunks = ChunkList.convertChromSizes(chrs, splitsize = splitsize, offset = args.flank)
        sets = chunks.split(items = bases/splitsize)
    else:
        chunks = ChunkList.read(args.bed, chromDict = chrs, min_offset = args.flank)
        sets = chunks.split(bases = bases)
    params = _PWMParameters(bam = args.bam, up = args.flank, down = args.flank, fasta = args.fasta,
                            lower = args.lower, upper = args.upper, atac = args.atac, sym = args.sym)
    pool = Pool(processes = args.cores)
    tmp = pool.map(_pwmHelper, zip(sets,itertools.repeat(params)))
    pool.close()
    pool.join()
    n = 0.0
    result = np.zeros((len(params.nucleotides), params.up + params.down + 1))
    for i in tmp:
        result += i[0]
        n += i[1]
    result /= n
    if args.bed:
        normfreqs = seq.getNucFreqsFromChunkList(chunks, args.fasta, params.nucleotides)
    else:
        normfreqs = seq.getNucFreqs(args.fasta, params.nucleotides)
    result = result / np.reshape(np.repeat(normfreqs,result.shape[1]),result.shape)
    if args.sym:
        #Symmetrize
        left = result[:,0:(args.flank + 1)]
        right = result[:,args.flank:]
        rightflipped = np.fliplr(np.flipud(right))
        combined = (left + rightflipped) / 2
        result = np.hstack((combined, np.fliplr(np.flipud(combined[:,0:args.flank]))))
    #save
    pwm = PWM(result, args.flank, args.flank, params.nucleotides)
    pwm.save(args.out + '.PWM.txt')

