#!/usr/bin/env python
"""
Classes for computing nucleosome occupancy

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""


from scipy import signal, optimize, stats
import numpy as np
import pysam
import matplotlib.pyplot as plt
import PyATAC as PA

class InsertDistribution:
    """Class for modelling insert size distribution"""
    def __init__(self,  i_lower = 0, i_upper =2000):
        self.i_lower = i_lower
        self.i_upper = i_upper
    def getInsertSizes(self,bam, bed_list, maxreads = 500000):
        """Use bamfile to get insertsize distribution

        Uses maxreads (or total reads if less) reads to determine
        """
        bamfile=pysam.Samfile(bam, "rb")
        i=0#counter
        sizes = np.zeros(self.i_upper - self.i_lower)
        # loop through reads
        for chunk in bed_list:
            for p2_rds in bamfile.fetch(chunk.chrom,chunk.start,chunk.end):
                #check mapping quality
                if p2_rds.mapq>=30 and p2_rds.is_proper_pair and not p2_rds.is_reverse:
                    # get insert size
                    ilen = abs(p2_rds.tlen)-8
                    if ilen<self.i_upper and ilen>=self.i_lower:
                        sizes[ilen - self.i_lower]+=1
                        #increase counter by 1
                        i+=1
            if i>maxreads:
                break
        bamfile.close()
        self.insertsizes = PA.InsertSizes(sizes / sum(sizes), self.i_lower, self.i_upper)
    def modelNFR(self, boundary = 115):
        """Model NFR distribution with exponential distribution"""
        b = np.where(self.insertsizes.get(self.i_lower,boundary) == max(self.insertsizes.get()))[0][0]+10 + self.i_lower
        def exp_pdf(x,*p): #defines the PDF
            k=p[0]
            a=p[1]
            x=x-b
            return a*k*np.exp(-k*x)
        x = np.array(range(b,boundary))
        p0 = (.1,1)
        coeff, var_matrix = optimize.curve_fit(exp_pdf,x, self.insertsizes.get(b,boundary),
                                               p0=p0)
        nfr = np.concatenate((self.insertsizes.get(self.i_lower,boundary), exp_pdf(np.array(range(boundary,self.i_upper)),*coeff)))
        nfr[nfr==0] = min(nfr[nfr!=0])*0.01
        self.nfr_fit = PA.InsertSizes(nfr,self.i_lower,self.i_upper)
        nuc = np.concatenate((np.zeros(boundary-self.i_lower),
                            self.insertsizes.get(boundary,self.i_upper) -
                            self.nfr_fit.get(boundary,self.i_upper)))
        nuc[nuc<=0]=min(min(nfr)*0.1,min(nuc[nuc>0])*0.001)
        self.nuc_fit = PA.InsertSizes(nuc, self.i_lower, self.i_upper)
    def plotFits(self,filename=None):
        """plot the Fits"""
        fig = plt.figure()
        plt.plot(range(self.i_lower,self.i_upper),self.insertsizes.get(),
                 label = "Observed")
        #plt.plot(range(self.i_lower,self.i_upper),self.smoothed.get(), label = "Smoothed")
        plt.plot(range(self.i_lower,self.i_upper),self.nuc_fit.get(), label = "Nucleosome Fit")
        plt.plot(range(self.i_lower,self.i_upper),self.nfr_fit.get(), label = "NFR Fit")
        plt.legend()
        plt.xlabel("Distance between Insertions (inclusive)")
        plt.ylabel("Relative Frequency")
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            #Also save text output!
            filename2 = ".".join(filename.split(".")[:-1]+['txt'])
            out = np.vstack((self.insertsizes.get(), #self.smoothed.get(),
                            self.nuc_fit.get(), self.nfr_fit.get()))
            np.savetxt(filename2,out,delimiter="\t")
        else:
            fig.show()

class OccupancyCalcParams:
    """Class with parameters for occupancy determination"""
    def __init__(self, i_lower, i_upper , insert_dist):
        self.i_lower = i_lower
        self.i_upper = i_upper
        self.smooth_mat = np.tile(signal.gaussian(151,25),(i_upper-i_lower,1))
        nuc_probs = insert_dist.nuc_fit.get(i_lower,i_upper)
        self.nuc_probs = nuc_probs /np.sum(nuc_probs)
        nfr_probs = insert_dist.nfr_fit.get(i_lower,i_upper)
        self.nfr_probs = nfr_probs /np.sum(nfr_probs)
        self.alphas = np.linspace(0, 1, 101)
        self.x = map(lambda alpha: np.log(alpha * self.nuc_probs + (1 - alpha) * self.nfr_probs), self.alphas)
        self.w = 75
        self.l = len(self.x)


def calculateOccupancy(inserts, params):
    """function to calculate occupancy based on insert distribution
    also takes OccupancyCalcParams as input
    """
    logliks = np.array(map(lambda j: np.sum(params.x[j]*inserts),range(params.l)))
    logliks[np.isnan(logliks)] = -float('inf')
    occ = params.alphas[np.argmax(logliks)]
    return occ



class OccupancyTrack(PA.Track):
    """Class for computing nucleosome occupancy"""
    def __init__(self, chrom, start, end):
        PA.Track.__init__(self, chrom, start, end, "occupancy")
    def calculateOccupancyMLE_smooth(self, mat, params):
        """Calculate Occupancy track"""
        offset=self.start - mat.start
        if offset<params.w:
            raise Exception("For calculateOccupancyMLE_smooth, mat does not have sufficient flanking regions"),offset
        self.vals=np.ones(self.end - self.start)*float('nan')
        for i in range(len(self.vals)):
            new_inserts = np.sum(mat.get(i_lower = params.i_lower,i_upper = params.i_upper,
                                         start = self.start+i-params.w, end = self.start+i+params.w+1)*params.smooth_mat,
                                         axis = 1)
            if sum(new_inserts)>0:
                self.vals[i] = calculateOccupancy(new_inserts,params)

class OccPeak:
    def __init__(self, pos, chunk):
        """Class for storing occupancy peaks"""
        self.chrom = chunk.chrom
        self.pos = pos
        self.occ = chunk.occ.get(pos = pos)
        self.reads = chunk.cov.get(pos = pos)
    def writePeak(self, handle):
        """write bed line for peak"""
        handle.write("\t".join(map(str,[self.chrom,self.pos,self.pos+1,self.occ,self.reads]))+"\n")

class OccupancyParams:
    """Class for storing parmeers related to Occupancy determination"""
    def __init__(self, occ_calc_params, sep = 120, min_occ = 0.25, min_reads = 1, flank = 60, 
                write_peaks = True, out = None, bam = None):
        self.occ_calc_params = occ_calc_params
        self.sep = sep
        self.min_occ = min_occ
        self.min_reads = min_reads
        self.flank = flank
        self.write_peaks = write_peaks
        self.out = out
        self.bam = bam

class OccChunk(PA.ChromChunk):
    """Class for calculating occupancy and occupancy peaks
    """
    def __init__(self, chrom, start, end):
        PA.ChromChunk.__init__(self, chrom, start, end)
        self.peaks = {}
    def calculateOcc(self, bamfile, params):
        """calculate occupancy for chunk"""
        self.reads = PA.ReadList(self.chrom, self.start - 75, self.end + 75)
        self.reads.extract_reads(bamfile, i_upper = params.i_upper)
        self.mat = PA.ReadMat2D(self.chrom, self.start - 75, self.end + 75 , 0, params.i_upper)
        self.mat.makeReadMat(self.reads)
        self.occ = OccupancyTrack(self.chrom,self.start,self.end)
        self.occ.calculateOccupancyMLE_smooth(self.mat, params)
    def getCov(self, i_upper):
        """Get read coveraged (smoothed in same way as occupancy) for regions"""
        self.cov = PA.CoverageTrack(self.chrom, self.start, self.end)
        self.cov.calculateCoverageSmooth(self.mat, 0, i_upper, 151, sd = 25)
    def callPeaks(self, sep, min_occ, min_reads):
        """Call peaks of occupancy profile"""
        peaks = PA.call_peaks(self.occ.vals, sep = sep, min_signal = min_occ)
        for peak in peaks:
            tmp = OccPeak(peak + self.start, self)
            if tmp.reads > min_reads:
                self.peaks[peak] = tmp
    def getNucDist(self, flank, i_upper):
        """Get nucleosomal insert distribution"""
        nuc_dist = np.zeros(i_upper)
        for peak in self.peaks.keys():
            sub = self.mat.get(start = self.peaks[peak].pos-flank,end = self.peaks[peak].pos+1+flank)
            sub_sum = np.sum(sub,axis=1)
            sub_sum = sub_sum / float(sum(sub_sum))
            nuc_dist += sub_sum
        return(nuc_dist)
    def process(self, bamfile, params):
        """proces chunk -- calculat occupancy, get coverage, call peaks"""
        self.calculateOcc(bamfile, params.occ_calc_params)
        self.getCov(params.occ_calc_params.i_upper)
        self.callPeaks(params.sep, params.min_occ, params.min_reads)
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self, name)
 


class OccChunkSet(PA.ChunkSet):
    """Class for storing sets of OccChunks
    """
    def __init__(self, name, chunks):
        PA.ChunkSet.__init__(self, name, chunks)
    def process(self, params):
        """Process occupancy chunk set -- get occupancy track, call peaks, get nuc_dist"""
        bamfile = pysam.Samfile(params.bam, "rb")
        nuc_dist = np.zeros(params.occ_calc_params.i_upper)
        handles = {}
        outputs = ['occ']
        if params.write_peaks:
            outputs.extend(['occ_peaks'])
        for i in outputs:
            handles[i] = open('tmp_nucleoatac/' + self.name + '.' + params.out + '.' +
                            i + '.tmp.bed','w')
        for chunk in self.chunks:
            chunk.process(bamfile, params)
            chunk.occ.write_track(handles['occ'])
            if params.write_peaks:
                for key in sorted(chunk.peaks.keys()):
                    chunk.peaks[key].writePeak(handles['occ_peaks'])
            nuc_dist += chunk.getNucDist(params.flank, params.occ_calc_params.i_upper)
            chunk.removeData()
        for key in handles.keys():
            handles[key].close()
        bamfile.close()
        return nuc_dist

def makeChunkSetList(bed_list, number):
    """turn list of ChromChunks into a list of ChunkSet (sets of ChromChunks)"""
    l = len(bed_list)
    bases = 0.0
    for i in range(l):
        bases += bed_list[i].end - bed_list[i].start
    if number >= l:
        sets = []
        z = int(np.log10(l)+1)
        for j in range(l):
            sets.append(OccChunkSet(str(j).zfill(z),bed_list[j:(j+1)]))
    else:
        x = min(number*3,l)
        n = bases/x
        sets = []
        z = int(np.log10(x)+1)
        included = 0
        i = 0
        for j in range(x):
            a = i
            while included < (j+1)*n and i < l:
                included += bed_list[i].end - bed_list[i].start
                i += 1
            b = i
            sets.append(OccChunkSet(str(j).zfill(z),bed_list[a:b]))
    return(sets)


def read_bed(bedfile, genome_dict, max_offset = 0):
        """Make a dictionary of Chromosome insances from a bedfile"""
        infile = open(bedfile,"r")
        out = []
        while 1:
            in_line = infile.readline().rstrip('\n').split("\t")
            if not in_line[0]:
                break
            chrom = in_line[0]
            start = int(in_line[1])
            if start < max_offset:
                start = max_offset
            end = int(in_line[2])
            if end > (genome_dict[chrom] - max_offset):
                end = genome_dict[chrom] - max_offset
            out.extend([OccChunk(chrom, start, end)])
        infile.close()
        return out



