
import pysam
import numpy as np


class BedGraphFile:
    def __init__(self,bedgraph):
        self.tbx = pysam.Tabixfile(bedgraph)
    def read(self, chrom, start, end, empty = np.nan):
        out = np.ones(end-start)*empty
        if chrom in self.tbx.contigs:
            for row in self.tbx.fetch(chrom,start,end, parser=pysam.asTuple()):
                out[max(int(row[1])-start,0):min(int(row[2])-start,end-start)] = float(row[3])
        return out
    def close(self):
        self.tbx.close()


def tabix_bedgraph(bedgraph):
    pysam.tabix_compress(bedgraph,bedgraph+'.gz')
    pysam.tabix_index(bedgraph+'.gz', seq_col=0, start_col=1, end_col=2, zerobased=True)



