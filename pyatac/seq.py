"""
General tools for dealing with ATAC-Seq data using Python.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import string
import numpy as np
import pysam

def get_sequence(chunk, fastafile):
    """obtain sequence for an interval

        chunk:  chunk object for which sequenceuence is to be fetched
        fastafile: filename for fasta file with sequenceuence
    """
    handle = pysam.FastaFile(fastafile)
    sequence = handle.fetch(chunk.chrom, chunk.start, chunk.end)
    if chunk.strand == "-":
        sequence = reverse_complement(sequence)
    handle.close()
    return sequence.upper()


DNA_Translation = string.maketrans('ACGT', 'TGCA')

def complement(sequence):
    """Get complement of DNA sequenceuence"""
    return sequence.translate(DNA_Translation)


def reverse_complement(sequence):
    """Get reverse complement of DNA sequenceuence"""
    return complement(sequence[::-1])


def seq_to_mat(sequence, nucleotides):
    """Turn sequenceuence into matrix encoding"""
    l = len(nucleotides[0])
    if not all(x==l for x in map(len,nucleotides)):
        raise Exception("Usage Error! Nucleotides must all be of same length! No mixing single nucleotides with dinucleotides, etc")
    mat = np.zeros((len(nucleotides),len(sequence)-l+1))
    for i in range(len(nucleotides)):
        mat[i] = np.array(map(int,[sequence[j:j+l] ==nucleotides[i] for j in range(len(sequence)-l+1)]))
    return(mat)

def getNucFreqs(fasta, nucleotides):
    """Get genomewide nucleotide frequencies"""
    out = np.zeros(len(nucleotides))
    n = 0.0
    f = open(fasta,'r')
    for line in f:
        if line[0]!='>':
            sequence = line.rstrip('\n').upper()
            out += [sequence.count(i) for i in nucleotides]
            n += len(sequence)
    f.close()
    return out/n


def getNucFreqsFromChunkList(chunks, fasta, nucleotides):
    """Get nucleotide frequences within regions of genome"""
    out = np.zeros(len(nucleotides))
    n = 0.0
    handle = pysam.FastaFile(fasta)
    for chunk in chunks:
        sequence = handle.fetch(chunk.chrom, chunk.start, chunk.end)
        sequence = sequence.upper()
        out += [sequence.count(i) for i in nucleotides]
        n += len(sequence)
    handle.close()
    return out/n

