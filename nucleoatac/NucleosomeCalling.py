"""
Script with classes and functions for nucleosome calling.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import numpy as np
from scipy import optimize, signal
from copy import copy
from bisect import bisect_left
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
from nucleoatac.multinomial_cov import calculateCov
from nucleoatac.Occupancy import OccupancyTrack
from pyatac.tracks import Track, CoverageTrack
from pyatac.chunk import Chunk
from pyatac.utils import call_peaks, reduce_peaks, read_chrom_sizes_from_bam
from pyatac.chunkmat2d import FragmentMat2D, BiasMat2D
from pyatac.bias import InsertionBiasTrack, PWM


#import warnings
#warnings.filterwarnings('error')


class SignalTrack(Track):
    """Class for getting V-plot signal"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "signal")
    def calculateSignal(self, mat, vmat):
        offset=self.start-mat.start-vmat.w
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate signal")
        self.vals = signal.correlate(mat.get(vmat.lower,vmat.upper,
                                              mat.start + offset, mat.end - offset),
                                       vmat.mat,mode = 'valid')[0]

class NormSignalTrack(Track):
    """Class for storing normalized signal track"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "normalized signal")
    def calculateNormSignal(self, raw, bias):
        self.vals = raw.get(self.start, self.end) - bias.get(self.start,self.end)

class BiasTrack(Track):
    """Class for getting Bias Signal Track-- Background model"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "bias")
    def calculateBackgroundSignal(self, mat, vmat, nuc_cov):
        offset=self.start-mat.start-vmat.w
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate signal")
        self.vmat = vmat
        self.bias_mat = mat
        self.cov = CoverageTrack(self.chrom, self.start, self.end)
        self.cov.calculateCoverage(self.bias_mat, vmat.lower,
                                   vmat.upper, vmat.w*2+1)
        self.nuc_cov = nuc_cov.vals
        self.vals = signal.correlate(self.bias_mat.get(vmat.lower,vmat.upper,
                                                         self.bias_mat.start + offset,
                                                         self.bias_mat.end - offset),
                                       vmat.mat,mode = 'valid')[0]
        self.vals = self.vals * self.nuc_cov/ self.cov.vals



class SignalDistribution:
    """Class for determining distribution of signal"""
    def __init__(self, position, vmat, bias_mat, reads):
        self.position = position
        self.reads = reads
        self.vmat = vmat
        bias_mat = bias_mat.get(vmat.lower,vmat.upper,position - vmat.w,position  + vmat.w + 1)
        self.prob_mat = bias_mat / np.sum(bias_mat)
        self.probs = self.prob_mat.flatten()
    def simulateReads(self):
        sim_vect = np.random.multinomial(self.reads,self.probs)
        sim_mat = np.reshape(sim_vect, self.vmat.mat.shape)
        return sim_mat
    def simulateDist(self, numiters = 1000):
        self.scores = map(lambda x: np.sum(self.simulateReads() * self.vmat.mat),range(numiters))
    def analStd(self):
        flatv = np.ravel(self.vmat.mat)
        var = calculateCov(self.probs, flatv, self.reads)
        return np.sqrt(var)
    def analMean(self):
        return np.sum(self.prob_mat * self.vmat.mat * self.reads)



def norm(x, v, w, mean):
    """compute values of normal pdf with given mean and sd at values in x"""
    norm = (1.0/(np.sqrt(2*np.pi*v)) *
        np.exp(-(x - mean)**2/(2*v)))
    norm = norm * (w/max(norm))
    return norm

class Nucleosome(Chunk):
    """Class for storing information about a single nucleosome"""
    def __init__(self, pos,nuctrack):
        self.chrom = nuctrack.chrom
        self.start = pos
        self.end = pos + 1
        self.nfr_cov = nuctrack.nfr_cov.get(pos = pos)
        self.nuc_cov = nuctrack.nuc_cov.get(pos = pos)
        self.nuc_signal = nuctrack.nuc_signal.get(pos = pos)
        self.norm_signal = nuctrack.norm_signal.get(pos = pos)
        self.smoothed = nuctrack.smoothed.get(pos= pos)
    def getLR(self,nuctrack):
        mat = nuctrack.mat.get(nuctrack.params.lower,nuctrack.params.upper,
                                self.start - nuctrack.params.vmat.w, self.start + nuctrack.params.vmat.w +1)
        null_mat = nuctrack.bias_mat.get(nuctrack.params.lower,nuctrack.params.upper,
                                self.start - nuctrack.params.vmat.w, self.start + nuctrack.params.vmat.w +1)
        bias_mat =nuctrack.bias_mat_prenorm.get(nuctrack.params.lower,nuctrack.params.upper,
                                self.start - nuctrack.params.vmat.w, self.start + nuctrack.params.vmat.w +1)
        nuc_model = nuctrack.params.vmat.mat * bias_mat
        nuc_model = nuc_model / np.sum(nuc_model)
        null_model = null_mat / np.sum(null_mat)
        nuc_lik = np.sum(np.log(nuc_model) * mat)
        null_lik = np.sum(np.log(null_model) * mat)
        self.lr = nuc_lik - null_lik
    def getZScore(self, nuctrack):
        s = SignalDistribution(self.start, nuctrack.params.vmat, nuctrack.bias_mat,
                                self.nuc_cov)
        std = s.analStd()
        self.z = self.norm_signal / std
    def getOcc(self, nuctrack):
        try:
            self.occ = nuctrack.occ.get(pos = self.start)
            self.occ_lower = nuctrack.occ_lower.get(pos = self.start)
            self.occ_upper = nuctrack.occ_upper.get(pos = self.start)
        except:
            self.occ = np.nan
            self.occ_lower = np.nan
            self.occ_upper = np.nan
    def getFuzz(self, nuctrack):
        def addNorms(x,params):
            """Add several normal distributions together"""
            l = len(x)
            fit = np.zeros(l)
            i = len(params)/3
            for j in range(i):
                fit += norm(x,params[j*3],params[3*j+1],params[3*j+2])
            return fit
        def err_func(pars,y):
            """error function for normal fit; to be used for fitNorm"""
            x = np.linspace(0,len(y)-1,len(y))
            return sum((addNorms(x, pars) - y)**2)
        def fitNorm(guess, bound, sig):
            """Fit a normal to the signal with lower and upperbounds to sd"""
            a = (sig,)
            res = optimize.minimize(err_func,guess,args = a, bounds=bound,method="L-BFGS-B")
            return res
        index = self.start - nuctrack.start
        allnucs = nuctrack.sorted_nuc_keys
        x = bisect_left(allnucs,index)
        if x == 0:
            left = index - nuctrack.params.nonredundant_sep/3
            means = (nuctrack.params.nonredundant_sep/3,)
        elif index - allnucs[x-1] < nuctrack.params.nonredundant_sep:
            left = allnucs[x-1]
            means = (index - allnucs[x-1],0)
        else:
            left = index - nuctrack.params.nonredundant_sep/3
            means = (nuctrack.params.nonredundant_sep/3,)
        if x == len(allnucs)-1:
            right = index + nuctrack.params.nonredundant_sep/3 + 1
        elif allnucs[x+1] - index < nuctrack.params.nonredundant_sep:
            right = allnucs[x+1]
            means += (allnucs[x+1] - left,)
        else:
            right = index + nuctrack.params.nonredundant_sep/3 +1
        sig = nuctrack.smoothed.vals[left:right]
        sig[sig<0] = 0
        if len(means)==1:
            bounds = ((2**2,50**2),(0.001,max(sig)*1.1),(means[0]-10,means[0]+10))
            guesses = (nuctrack.params.smooth_sd ** 2,max(sig)*0.9,means[0])
        elif len(means)==2:
            bounds = ((2**2,50**2),(0.001,max(sig)*1.1),(means[0]-10,means[0]+10),
                        (2**2,50**2),(0.001,max(sig)*1.1),(means[1]-10,means[1]+10))
            guesses = (nuctrack.params.smooth_sd ** 2,max(sig)*0.9,means[0],
                    nuctrack.params.smooth_sd ** 2,max(sig)*0.9,means[1])
        elif len(means)==3:
            bounds = ((2**2,50**2),(0.001,max(sig)*1.1),(means[0]-10,means[0]+10),
                        (2**2,50**2),(0.001,max(sig)*1.1),(means[1]-10,means[1]+10),
                        (2**2,50**2),(0.001,max(sig)*1.1),(means[2]-10,means[2]+10))
            guesses = (nuctrack.params.smooth_sd ** 2,max(sig)*0.9,means[0],
                    nuctrack.params.smooth_sd ** 2,max(sig)*0.9,means[1],
                    nuctrack.params.smooth_sd ** 2,max(sig)*0.9,means[2])
        res= fitNorm(guesses, bounds, sig)
        self.fuzz= np.sqrt(res['x'][0])
        self.weight = res['x'][1]
        self.fit_pos = res['x'][2]+left
    def asBed(self):
        out = "\t".join(map(str,[self.chrom, self.start, self.end, self.z, self.occ, self.occ_lower, self.occ_upper, self.lr,
               self.norm_signal, self.nuc_signal, self.nuc_cov, self.nfr_cov,
                self.fuzz]))
        return out
    def write(self, handle):
       handle.write(self.asBed() + "\n")


class NucParameters:
    """Class for storing parameters related to nucleosome calling"""
    def __init__(self, vmat, fragmentsizes, bam, fasta, pwm,
                 occ_track = None, atac = True,
                 sd = 25, nonredundant_sep = 120, redundant_sep = 25,
                 min_z = 3, min_lr = 0, min_reads = 1):
        self.atac = atac
        self.vmat = vmat
        self.lower = vmat.lower
        self.upper= vmat.upper
        self.window = vmat.mat.shape[1]
        self.fragmentsizes= fragmentsizes
        self.min_reads = min_reads
        self.min_z = min_z
        self.min_lr = min_lr
        self.smooth_sd = sd
        self.redundant_sep = redundant_sep
        self.nonredundant_sep = nonredundant_sep
        self.fasta = fasta
        self.pwm = PWM.open(pwm)
        self.chrs = read_chrom_sizes_from_bam(bam)
        self.bam = bam
        self.occ_track = occ_track



class NucChunk(Chunk):
    """Class for storing and determining collection of nucleosome positions
    """
    def __init__(self, chunk):
        self.start = chunk.start
        self.end = chunk.end
        self.chrom = chunk.chrom
    def initialize(self, parameters):
        self.params = parameters
    def getFragmentMat(self):
        self.mat = FragmentMat2D(self.chrom, self.start - max(self.params.window,self.params.upper/2+1),
                                 self.end + max(self.params.window,self.params.upper/2+1), 0, self.params.upper, atac = self.params.atac)
        self.mat.makeFragmentMat(self.params.bam)
    def makeBiasMat(self):
        self.bias_mat = BiasMat2D(self.chrom, self.start - self.params.window,
                                 self.end + self.params.window, 0, self.params.upper)
        bias_track = InsertionBiasTrack(self.chrom, self.start - self.params.window - self.params.upper/2,
                                  self.end + self.params.window + self.params.upper/2 + 1, log = True)
        if self.params.fasta is not None:
            bias_track.computeBias(self.params.fasta, self.params.chrs, self.params.pwm)
            self.bias_mat.makeBiasMat(bias_track)
        self.bias_mat_prenorm = BiasMat2D(self.chrom, self.start - self.params.window,
                                 self.end + self.params.window, 0, self.params.upper)
        self.bias_mat_prenorm.mat = copy(self.bias_mat.mat)
        self.bias_mat.normByInsertDist(self.params.fragmentsizes)
    def getNucSignal(self):
        """Gets Nucleosome Signal Track"""
        self.nuc_cov = CoverageTrack(self.chrom, self.start,
                                     self.end)
        self.nuc_cov.calculateCoverage(self.mat, self.params.lower, self.params.upper,
                                        self.params.window)
        self.bias = BiasTrack(self.chrom, self.start,
                                     self.end)
        self.bias.calculateBackgroundSignal(self.bias_mat, self.params.vmat, self.nuc_cov)
        self.nuc_signal = SignalTrack(self.chrom, self.start,
                                     self.end)
        self.nuc_signal.calculateSignal(self.mat, self.params.vmat)
        self.norm_signal = NormSignalTrack(self.chrom, self.start, self.end)
        self.norm_signal.calculateNormSignal(self.nuc_signal,self.bias)
    def getNFR(self):
        """get number of reads of sub-nucleosomal length"""
        self.nfr_cov = CoverageTrack(self.chrom, self.start, self.end)
        self.nfr_cov.calculateCoverage(self.mat, 0, self.params.lower,
                                                        self.params.window)
    def smoothSignal(self):
        """Smooth thenormalized signal track"""
        window_len = 6 * self.params.smooth_sd + 1
        self.smoothed = Track(self.chrom,self.start,self.end, "Smooth Signal")
        tmp = copy(self.norm_signal.vals)
        self.smoothed.assign_track(tmp)
        self.smoothed.vals[ self.smoothed.vals < 0] = 0
        self.smoothed.smooth_track(window_len, window = "gaussian",
                             sd = self.params.smooth_sd, mode = 'same',
                             norm = True)
    def getOcc(self):
        """gets occupancy track-- either reads in from bw handle given, or makes new"""
        self.occ = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ.read_track(self.params.occ_track)
        lower_file = self.params.occ_track[:-11] + 'lower_bound.bedgraph.gz'
        self.occ_lower = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ_lower.read_track(lower_file)
        upper_file = self.params.occ_track[:-11] + 'upper_bound.bedgraph.gz'
        self.occ_upper = Track(self.chrom,self.start,self.end,"Occupancy")
        self.occ_upper.read_track(upper_file)
    def findAllNucs(self):
        """Find peaks in data"""
        self.nuc_collection = {}
        combined = self.norm_signal.vals + self.smoothed.vals
        #find peaks in normalized sigal
        cands1 = call_peaks(combined, min_signal = 0,
                                sep = self.params.redundant_sep,
                                boundary = self.params.nonredundant_sep/2, order = self.params.redundant_sep/2)
        for i in cands1:
            nuc = Nucleosome(i + self.start, self)
            if nuc.nuc_cov > self.params.min_reads:
                nuc.getLR(self)
                if nuc.lr > self.params.min_lr:
                    nuc.getZScore(self)
                    if nuc.z >= self.params.min_z:
                        nuc.getOcc(self)
                        self.nuc_collection[i] = nuc
        self.sorted_nuc_keys = np.array(sorted(self.nuc_collection.keys()))
        self.nonredundant = reduce_peaks( self.sorted_nuc_keys,
                                            map(lambda x: self.nuc_collection[x].z, self.sorted_nuc_keys),
                                                self.params.nonredundant_sep)
        self.redundant = np.setdiff1d(self.sorted_nuc_keys, self.nonredundant)
    def fit(self):
        x = np.linspace(0,self.length() -1, self.length())
        fit = np.zeros(self.length())
        for nuc in self.sorted_nuc_keys:
            self.nuc_collection[nuc].getFuzz(self)
            fit += norm(x,self.nuc_collection[nuc].fuzz**2, self.nuc_collection[nuc].weight, self.nuc_collection[nuc].fit_pos)
        self.fitted = Track(self.chrom, self.start, self.end,
                            "Fitted Nucleosome Signal")
        self.fitted.assign_track(fit)
    def makeInsertionTrack(self):
        """make insertion track for chunk"""
        self.ins = self.mat.getIns()
    def process(self, params):
        """wrapper to carry out all methods needed to call nucleosomes and nfrs"""
        self.initialize(params)
        self.getFragmentMat()
        self.makeBiasMat()
        self.getNucSignal()
        self.getNFR()
        self.smoothSignal()
        if params.occ_track is not None:
            self.getOcc()
        self.findAllNucs()
        self.fit()
        self.makeInsertionTrack()
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self,name)

