#!/usr/bin/env python
"""
General tools for dealing with ATAC-Seq data using Python.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import gzip
import warnings

class Chunk():
    """Class that stores reads for a particular chunk of the genome"""
    def __init__(self, chrom, start, end, weight = 1, name = "region", strand = "*"):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.weight = weight
        self.strand = strand
        self.name = name
    def length(self):
        return self.end - self.start
    def asBed(self):
        """represent output as bed"""
        out = "\t".join(map(str,[self.chrom,self.start,self.end,self.weight,self.name,self.strand]))
        return out
    def slop(self, chromDict, up = 0, down = 0, new = False):
        """extend region, checking for chromosomal constraints"""
        if self.strand == "-":
            newStart = max(0,self.start - down)
            newEnd = min(chromDict[self.chrom],self.end + up)
        else:
            newStart = max(0,self.start - up)
            newEnd = min(chromDict[self.chrom],self.end + down)
        if new:
            out = Chunk(self.chrom, newStart, newEnd,
                        weight = self.weight, name = self.name, strand = self.strand)
            return out
        else:
            self.start = newStart
            self.end = newEnd
    def center(self, new = False):
        if self.strand == "-":
            newEnd = self.end - (self.length()/2)
            newStart = newEnd - 1
        else:
            newStart = self.start + (self.length()/2)
            newEnd = newStart +1
        if new:
            out = Chunk(self.chrom, newStart, newEnd,
                        weight = self.weight, name = self.name, strand = self.strand)
            return out
        else:
            self.start = newStart
            self.end = newEnd


def _chunkCompare(chunk1, chunk2):
    """Compare positions of two chunks"""
    if chunk1.chrom < chunk2.chrom:
        return -1
    elif chunk1.chrom > chunk2.chrom:
        return 1
    else:
        if chunk1.start < chunk2.start:
            return -1
        elif chunk2.start > chunk2.start:
            return 1
        else:
            return 0

class ChunkList(list):
    def __init__(self, *args):
        list.__init__(self, args)
    def extend(self, *args):
        if len(args)!=1:
            raise ValueError("Wrong number of arguments")
        elif isinstance(args[0], ChunkList):
            list.extend(self, args[0])
        else:
            raise ValueError("Expecting ChunkList")
    def append(self, *args):
        if len(args)!=1:
            raise ValueError("Wrong number of arguments")
        elif isinstance(args[0],Chunk):
            list.append(self, args[0])
        else:
            raise ValueError("Expecting Chunk")
    def insert(self, *args):
        if len(args)!=1:
            raise ValueError("Wrong number of arguments")
        elif isinstance(args[1],Chunk):
            list.insert(self, args[0], args[1])
        else:
            raise ValueError("Expecting Chunk")
    def sort(self):
        """sort regions"""
        list.sort(self, cmp = _chunkCompare)
    def isSorted(self):
        """check that regions are sorted"""
        return all([_chunkCompare(self[i],self[i+1])==-1 for i in xrange(len(self)-1)])
    def slop(self, chromDict, up = 0, down = 0, new = False):
        out = ChunkList()
        for i in self:
            out.append(i.slop(chromDict, up, down, new = True))
        if new:
            return out
        else:
            self[:] = out
    def merge(self, new = False, sep = -1):
        """Merge overlapping or nearby regions"""
        if not self.isSorted():
            self.sort()
        out = ChunkList()
        previous = self[0]
        for i in range(1,len(self)):
            if self[i].chrom == previous.chrom and self[i].start <= (previous.end + sep):
                previous.end = max(self[i].end, previous.end)
            else:
                out.append(previous)
                previous = self[i]
        out.append(previous)
        if new:
            return out
        else:
            self[:] = out
    def asBed(self):
        """format regions as bed"""
        out = ""
        for chunk in self:
            out += chunk.asBed() +"\n"
        return out
    @staticmethod
    def read(bedfile, weight_col=None, strand_col = None, name_col = None, chromDict = None,
                min_offset = None, min_length = 1, chrom_source = "FASTA file"):
        """Make a list of chunks from a tab-delimited bedfile"""
        if bedfile[-3:] == '.gz':
            infile = gzip.open(bedfile,"r")
        else:
            infile = open(bedfile,"r")
        out = ChunkList()
        weight = None
        strand = "+"
        name = None
        if chromDict is not None:
            bad_chroms = []
        for line in infile:
            in_line = line.rstrip('\n').split("\t")
            if weight_col:
                weight=in_line[weight_col-1]
            if strand_col:
                strand = in_line[strand_col-1]
            if name_col:
                name = in_line[strand_col-1]
            start = int(in_line[1])
            end = int(in_line[2])
            chrom = in_line[0]
            if chromDict is not None and chrom not in chromDict.keys():
                bad_chroms.append(chrom)
                continue
            if min_offset:
                if start < min_offset:
                    start = min_offset
                if end > (chromDict[in_line[0]] - min_offset):
                    end = chromDict[in_line[0]] - min_offset
            if end - start >= min_length:
                out.append(Chunk(in_line[0],start, end,
                                   weight = weight, strand = strand, name = name))
        infile.close()
        if chromDict is not None and len(bad_chroms)>0:
            bad_chroms = set(bad_chroms)
            warn_message = (str(len(bad_chroms)) + " chromosome names in bed file not included in " +
                            chrom_source + ":\n" + "\n".join(bad_chroms) +  "\n " + 
                            "These regions will be ignored in subsequent analysis")
            warnings.warn(warn_message)
        return out
    @staticmethod
    def convertChromSizes(chromDict, splitsize = None, offset = 0):
        """Convert dictionary of chromosome sizes to Chunklist"""
        chrs = sorted(chromDict.keys())
        if splitsize is None:
            return ChunkList(*(Chunk(chrom, offset, chromDict[chrom]-offset) for chrom in chrs))
        else:
            out = ChunkList()
            for chrom in chrs:
                out.extend(ChunkList(*(Chunk(chrom, i, min(i + splitsize, chromDict[chrom] - offset))
                        for i in xrange(offset, chromDict[chrom] - offset, splitsize))))
            return out
    def split(self, bases = None, items = None):
        """splits list of chunks into set of sublists"""
        if bases is not None:
            out = []
            i = 0
            j = 0
            k = 0
            for k in range(len(self)):
                j += self[k].length()
                if j > bases:
                    out.append(self[i:(k+1)])
                    j = 0
                    i = k+1
            if k >= i:
                out.append(self[i:(k+1)])
            return out
        elif items is not None:
            out = [ self[i:i+items] for i in xrange(0,len(self),items)]
            return out
        else:
            raise Exception("Need to provide items or bases argument!")
    def checkChroms(self, chroms, chunklist_source = "bed file", chrom_source = "fasta file", warn = "Regions on these chromosomes will be ignored in analysis"):
        bad_chroms = set([x.chrom for x in self if x.chrom not in chroms])
        if len(bad_chroms)>0:
            self[:] = [x for x in self if x.chrom in chroms]
            warn_message = (str(len(bad_chroms)) + " chromosome names in " + chunklist_source + " not included in "
                            + chrom_source + ":\n" + "\n".join(bad_chroms) +  "\n " + warn)
            warnings.warn(warn_message)


