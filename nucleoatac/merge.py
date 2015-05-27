"""
Script to merge nuc positions

@author: Alicia Schep
"""
##### IMPORT MODULES #####
#

from pyatac.utils import shell_command
from pyatac.chunk import Chunk, ChunkList
import gzip
import pysam
import os

class MergedNuc(Chunk):
    def __init__(self, chrom, start, end, occ, occ_lower, occ_upper, reads, source):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.occ = occ
        self.occ_lower = occ_lower
        self.occ_upper = occ_upper
        self.reads = reads
        self.source = source
    def asBed(self):
        out = "\t".join(map(str,[self.chrom,self.start,self.end,self.occ,self.occ_lower,self.occ_upper,self.reads, self.source])) 
        return out
    def write(self, handle):
        """write bed line for peak"""
        handle.write(self.asBed() + "\n")

class NucList(ChunkList):
    def __init__(self, *args):
        list.__init__(self,args)
    @staticmethod
    def read(bedfile, source, min_occ = 0):
        """Make a list of chunks from a tab-delimited bedfile"""
        if bedfile[-3:] == '.gz':
            infile = gzip.open(bedfile,"r")
        else:
            infile = open(bedfile,"r")
        out = NucList()
        for line in infile:
            in_line = line.rstrip('\n').split("\t")
            start = int(in_line[1])
            end = int(in_line[2])
            if source == "occ":
                occ = float(in_line[3])
                occ_lower = float(in_line[4])
                occ_upper = float(in_line[5])
                reads = float(in_line[6])
            elif source == "nuc":
                occ = float(in_line[4])
                occ_lower = float(in_line[5])
                occ_upper = float(in_line[6])
                reads = float(in_line[10]) + float(in_line[11])
            else:
                raise Exception("source must be 'occ' or 'nuc'")
            if occ_lower >= min_occ:
                out.append(MergedNuc(in_line[0],start, end, occ, occ_lower, occ_upper, reads, source))
        infile.close()
        return out




def merge(occ_peaks, nuc_calls, sep = 120):
    keep = NucList()
    i = 0 #index for occ peaks
    j = 0 #index for nuc calls
    while i < len(occ_peaks) and j < len(nuc_calls):
        if occ_peaks[i].chrom < nuc_calls[j].chrom:
            keep.append(occ_peaks[i])
            i += 1
        elif occ_peaks[i].chrom > nuc_calls[j].chrom:
            keep.append(nuc_calls[j])
            j += 1
        elif occ_peaks[i].start < (nuc_calls[j].start - sep):
            keep.append(occ_peaks[i])
            i += 1
        elif occ_peaks[i].start > (nuc_calls[j].start + sep):
            keep.append(nuc_calls[j])
            j += 1
        else:
            i += 1
    while j < len(nuc_calls):
        keep.append(nuc_calls[j])
        j += 1
    while i < len(occ_peaks):
        keep.append(occ_peaks[i])
        i += 1
    return keep

def run_merge(args):
    if not args.out:
        args.out = '.'.join(os.path.basename(args.nucpos).split('.')[0:-3])
    occ = NucList.read(args.occpeaks, "occ", args.min_occ)
    nuc = NucList.read(args.nucpos, "nuc", args.min_occ)
    new = merge(occ, nuc, args.sep)
    out = open(args.out + '.nucmap_combined.bed','w')
    out.write(new.asBed())
    out.close()
    pysam.tabix_compress(args.out + '.nucmap_combined.bed', args.out + '.nucmap_combined.bed.gz',force = True)
    shell_command('rm ' + args.out + '.nucmap_combined.bed')
    pysam.tabix_index(args.out + '.nucmap_combined.bed.gz', preset = "bed", force = True)
 



