"""
Classes for computing nucleosome occupancy

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""


from scipy import signal, optimize, stats
import numpy as np
import matplotlib.pyplot as plt
import pyximport; pyximport.install()
from pyatac.fragmentsizes import FragmentSizes
from pyatac.tracks import Track, CoverageTrack
from pyatac.chunk import Chunk
from pyatac.utils import call_peaks, read_chrom_sizes_from_fasta
from pyatac.chunkmat2d import FragmentMat2D, BiasMat2D
from pyatac.bias import InsertionBiasTrack, PWM


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
    def __init__(self, lower, upper , insert_dist, ci = 0.95):
        self.lower = lower
        self.upper = upper
        #self.smooth_mat = np.tile(signal.gaussian(151,25),(upper-lower,1))
        nuc_probs = insert_dist.nuc_fit.get(lower,upper)
        self.nuc_probs = nuc_probs /np.sum(nuc_probs)
        nfr_probs = insert_dist.nfr_fit.get(lower,upper)
        self.nfr_probs = nfr_probs /np.sum(nfr_probs)
        self.alphas = np.linspace(0, 1, 101)
        #self.x = map(lambda alpha: np.log(alpha * self.nuc_probs + (1 - alpha) * self.nfr_probs), self.alphas)
        self.l = len(self.alphas)
        self.cutoff = stats.chi2.ppf(ci,1)

def calculateOccupancy(inserts, bias, params):
    """function to calculate occupancy based on insert distribution
    also takes OccupancyCalcParams as input
    """
    nuc_probs = params.nuc_probs * bias
    nuc_probs = nuc_probs / np.sum(nuc_probs)
    nfr_probs = params.nfr_probs * bias
    nfr_probs = nfr_probs / np.sum(nfr_probs)
    x = map(lambda alpha: np.log(alpha * nuc_probs + (1 - alpha) * nfr_probs), params.alphas)
    logliks = np.array(map(lambda j: np.sum(x[j]*inserts),range(params.l)))
    logliks[np.isnan(logliks)] = -float('inf')
    occ = params.alphas[np.argmax(logliks)]
    #Compute upper and lower bounds for 95% confidence interval
    ratios = 2*(max(logliks)-logliks)
    lower = params.alphas[min(np.where(ratios < params.cutoff)[0])]
    upper = params.alphas[max(np.where(ratios < params.cutoff)[0])]
    return occ, lower, upper



class OccupancyTrack(Track):
    """Class for computing nucleosome occupancy"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "occupancy")
    def calculateOccupancyMLE(self, mat, bias_mat, params):
        """Calculate Occupancy track"""
        offset=self.start - mat.start
        if offset<params.flank:
            raise Exception("For calculateOccupancyMLE, mat does not have sufficient flanking regions"),offset
        self.vals=np.ones(self.end - self.start)*float('nan')
        self.lower_bound = np.ones(self.end - self.start)*float('nan')
        self.upper_bound =np.ones(self.end - self.start)*float('nan')
        for i in range(len(self.vals)):
            new_inserts = np.sum(mat.get(lower = 0, upper = params.upper,
                                         start = self.start+i-params.flank, end = self.start+i+params.flank+1),
                                         axis = 1)
            new_bias = np.sum(bias_mat.get(lower = 0, upper = params.upper,
                                         start = self.start+i-params.flank, end = self.start+i+params.flank+1),
                                         axis = 1)
            if sum(new_inserts)>0:
                self.vals[i],self.lower_bound[i],self.upper_bound[i] = calculateOccupancy(new_inserts, new_bias, params.occ_calc_params)

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
    def __init__(self, insert_dist, upper, fasta, pwm, sep = 120, min_occ = 0.25, min_reads = 1, flank = 30,
                write_peaks = True, out = None, bam = None):
        self.sep = sep
        self.chrs = read_chrom_sizes_from_fasta(fasta)
        self.fasta = fasta
        self.pwm = PWM.open(pwm)
        self.window = flank * 2 + 1
        self.min_occ = min_occ
        self.min_reads = min_reads
        self.flank = flank
        self.write_peaks = write_peaks
        self.out = out
        self.bam = bam
        self.upper = upper
        self.occ_calc_params = OccupancyCalcParams(0, upper, insert_dist)

class OccChunk(Chunk):
    """Class for calculating occupancy and occupancy peaks
    """
    def __init__(self, chunk):
        self.start = chunk.start
        self.end = chunk.end
        self.chrom = chunk.chrom
        self.peaks = {}
    def getFragmentMat(self):
        self.mat = FragmentMat2D(self.chrom, self.start - self.params.flank,
                                 self.end + self.params.flank, 0, self.params.upper)
        self.mat.makeFragmentMat(self.params.bam)
    def makeBiasMat(self):
        self.bias_mat = BiasMat2D(self.chrom, self.start - self.params.flank,
                                 self.end + self.params.flank, 0, self.params.upper)
        bias_track = InsertionBiasTrack(self.chrom, self.start - self.params.window - self.params.upper/2,
                                  self.end + self.params.window + self.params.upper/2 + 1, log = True)
        bias_track.computeBias(self.params.fasta, self.params.chrs, self.params.pwm)
        self.bias_mat.makeBiasMat(bias_track)
    def calculateOcc(self):
        """calculate occupancy for chunk"""
        self.occ = OccupancyTrack(self.chrom,self.start,self.end)
        self.occ.calculateOccupancyMLE(self.mat, self.bias_mat, self.params)
    def getCov(self):
        """Get read coverage for regions"""
        self.cov = CoverageTrack(self.chrom, self.start, self.end)
        self.cov.calculateCoverage(self.mat, 0, self.params.upper, self.params.window)
    def callPeaks(self):
        """Call peaks of occupancy profile"""
        peaks = call_peaks(self.occ.vals, sep = self.params.sep, min_signal = self.params.min_occ)
        for peak in peaks:
            tmp = OccPeak(peak + self.start, self)
            if tmp.reads > self.params.min_reads:
                self.peaks[peak] = tmp
    def getNucDist(self):
        """Get nucleosomal insert distribution"""
        nuc_dist = np.zeros(self.params.upper)
        for peak in self.peaks.keys():
            sub = self.mat.get(start = self.peaks[peak].pos-self.params.flank, end = self.peaks[peak].pos+1+self.params.flank)
            sub_sum = np.sum(sub,axis=1)
            sub_sum = sub_sum / float(sum(sub_sum))
            nuc_dist += sub_sum
        return(nuc_dist)
    def process(self, params):
        """proces chunk -- calculat occupancy, get coverage, call peaks"""
        self.params = params
        self.getFragmentMat()
        self.makeBiasMat()        
        self.calculateOcc()
        self.getCov()
        self.callPeaks()
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self, name)




