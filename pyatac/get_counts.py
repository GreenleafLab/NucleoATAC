"""
Script to get read counts  distribution

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary for python
import os
from pyatac.chunk import ChunkList
from pysam import AlignmentFile
import numpy as np

def _between(x,start,end):
    if x >= start and x < end:
        return True
    else:
        return False

def get_counts(args):
    """function to get fragment sizes

    """
    if args.out is None:
        args.out = '.'.join(os.path.basename(args.bed).split('.')[0:-1])  
    chunks = ChunkList.read(args.bed)
    mat = np.zeros(len(chunks), dtype=np.int)
    bamHandle = AlignmentFile(args.bam)
    j = 0
    for chunk in chunks:
        for read in bamHandle.fetch(chunk.chrom, max(0, chunk.start - args.upper), chunk.end + args.upper):
            if read.is_proper_pair and not read.is_reverse:
                if args.atac:
                    #get left position
                    l_pos = read.pos + 4
                    #get insert size
                    #correct by 8 base pairs to be inserion to insertion
                    ilen = abs(read.template_length) - 8
                else:
                    l_pos = read.pos
                    ilen = abs(read.template_length)
                r_pos = l_pos + ilen - 1
                if _between(ilen, args.lower, args.upper) and (_between(l_pos, chunk.start, chunk.end) or _between(r_pos, chunk.start, chunk.end)):
                    mat[j] += 1
        j += 1
    bamHandle.close()
    np.savetxt(args.out + ".counts.txt.gz", mat, delimiter="\n", fmt='%i')



