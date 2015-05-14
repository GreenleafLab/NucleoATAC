"""
Classes for computing nucleosome occupancy

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""


from scipy import signal, optimize
import numpy as np
import matplotlib.pyplot as plt
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
from pyatac.fragmentsizes import FragmentSizes
from pyatac.tracks import Track, CoverageTrack
from pyatac.chunk import Chunk
from pyatac.utils import call_peaks
from pyatac.chunkmat2d import FragmentMat2D

class FragmentMixDistribution:
    """Class for modelling insert size distribution"""
    def __init__(self,  lower = 0, upper =2000):
        self.lower = lower
        self.upper = upper
    def getFragmentSizes(self, bamfile, chunklist = None):
        self.fragmentsizes = FragmentSizes(self.lower, self.upper)
        self.fragmentsizes.calculateSizes(bamfile, chunks = chunklist)
    def modelNFR(self, boundary = 115):
        """Model NFR distribution with exponential distribution"""
        b = np.where(self.fragmentsizes.get(self.lower,boundary) == max(self.fragmentsizes.get(self.lower,boundary)))[0][0]+10 + self.lower
        def exp_pdf(x,*p): #defines the PDF
            k=p[0]
            a=p[1]
            x=x-b
            return a*k*np.exp(-k*x)
        x = np.array(range(b,boundary))
        p0 = (.1,1)
        coeff, var_matrix = optimize.curve_fit(exp_pdf,x, self.fragmentsizes.get(b,boundary),
                                               p0=p0)
        nfr = np.concatenate((self.fragmentsizes.get(self.lower,boundary), exp_pdf(np.array(range(boundary,self.upper)),*coeff)))
        nfr[nfr==0] = min(nfr[nfr!=0])*0.01
        self.nfr_fit = FragmentSizes(self.lower,self.upper, vals = nfr)
        nuc = np.concatenate((np.zeros(boundary-self.lower),
                            self.fragmentsizes.get(boundary,self.upper) -
                            self.nfr_fit.get(boundary,self.upper)))
        nuc[nuc<=0]=min(min(nfr)*0.1,min(nuc[nuc>0])*0.001)
        self.nuc_fit = FragmentSizes(self.lower, self.upper, vals = nuc)
    def plotFits(self,filename=None):
        """plot the Fits"""
        fig = plt.figure()
        plt.plot(range(self.lower,self.upper),self.fragmentsizes.get(),
                 label = "Observed")
        #plt.plot(range(self.lower,self.upper),self.smoothed.get(), label = "Smoothed")
        plt.plot(range(self.lower,self.upper),self.nuc_fit.get(), label = "Nucleosome Fit")
        plt.plot(range(self.lower,self.upper),self.nfr_fit.get(), label = "NFR Fit")
        plt.legend()
        plt.xlabel("Fragment size")
        plt.ylabel("Relative Frequency")
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            #Also save text output!
            filename2 = ".".join(filename.split(".")[:-1]+['txt'])
            out = np.vstack((self.fragmentsizes.get(), #self.smoothed.get(),
                            self.nuc_fit.get(), self.nfr_fit.get()))
            np.savetxt(filename2,out,delimiter="\t")
        else:
            fig.show()

class OccupancyCalcParams:
    """Class with parameters for occupancy determination"""
    def __init__(self, lower, upper , insert_dist):
        self.lower = lower
        self.upper = upper
        self.smooth_mat = np.tile(signal.gaussian(151,25),(upper-lower,1))
        nuc_probs = insert_dist.nuc_fit.get(lower,upper)
        self.nuc_probs = nuc_probs /np.sum(nuc_probs)
        nfr_probs = insert_dist.nfr_fit.get(lower,upper)
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



class OccupancyTrack(Track):
    """Class for computing nucleosome occupancy"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "occupancy")
    def calculateOccupancyMLE_smooth(self, mat, params):
        """Calculate Occupancy track"""
        offset=self.start - mat.start
        if offset<params.w:
            raise Exception("For calculateOccupancyMLE_smooth, mat does not have sufficient flanking regions"),offset
        self.vals=np.ones(self.end - self.start)*float('nan')
        for i in range(len(self.vals)):
            new_inserts = np.sum(mat.get(lower = params.lower,upper = params.upper,
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
    def write(self, handle):
        """write bed line for peak"""
        handle.write("\t".join(map(str,[self.chrom,self.pos,self.pos+1,self.occ,self.reads]))+"\n")


class OccupancyParameters:
    """Class for storing parmeers related to Occupancy determination"""
    def __init__(self, bed_list, insert_dist, upper, sep = 120, min_occ = 0.25, min_reads = 1, flank = 60,
                write_peaks = True, out = None, bam = None):
        self.sep = sep
        self.min_occ = min_occ
        self.min_reads = min_reads
        self.flank = flank
        self.write_peaks = write_peaks
        self.out = out
        self.bam = bam
        self.occ_calc_params = OccupancyCalcParams(0, upper, insert_dist)

class OccChunk(Chunk):
    """Class for calculating occupancy and occupancy peaks
    """
    def __init__(self, chunk):
        self.start = chunk.start
        self.end = chunk.end
        self.chrom = chunk.chrom
        self.peaks = {}
    def calculateOcc(self, bamfile, params):
        """calculate occupancy for chunk"""
        self.mat = FragmentMat2D(self.chrom, self.start - 75, self.end + 75 , 0, params.upper)
        self.mat.makeFragmentMat(bamfile)
        self.occ = OccupancyTrack(self.chrom,self.start,self.end)
        self.occ.calculateOccupancyMLE_smooth(self.mat, params)
    def getCov(self, upper):
        """Get read coveraged (smoothed in same way as occupancy) for regions"""
        self.cov = CoverageTrack(self.chrom, self.start, self.end)
        self.cov.calculateCoverageSmooth(self.mat, 0, upper, 151, sd = 25)
    def callPeaks(self, sep, min_occ, min_reads):
        """Call peaks of occupancy profile"""
        peaks = call_peaks(self.occ.vals, sep = sep, min_signal = min_occ)
        for peak in peaks:
            tmp = OccPeak(peak + self.start, self)
            if tmp.reads > min_reads:
                self.peaks[peak] = tmp
    def getNucDist(self, flank, upper):
        """Get nucleosomal insert distribution"""
        nuc_dist = np.zeros(upper)
        for peak in self.peaks.keys():
            sub = self.mat.get(start = self.peaks[peak].pos-flank,end = self.peaks[peak].pos+1+flank)
            sub_sum = np.sum(sub,axis=1)
            sub_sum = sub_sum / float(sum(sub_sum))
            nuc_dist += sub_sum
        return(nuc_dist)
    def process(self, params):
        """proces chunk -- calculat occupancy, get coverage, call peaks"""
        self.calculateOcc(params.bam, params.occ_calc_params)
        self.getCov(params.occ_calc_params.upper)
        self.callPeaks(params.sep, params.min_occ, params.min_reads)
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self, name)




