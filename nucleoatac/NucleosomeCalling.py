#!/usr/bin/env python
"""
Script with classes and functions for nucleosome calling.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import numpy as np
from scipy import optimize, signal
import PyATAC as PA
import Occupancy as Occ
from copy import copy
from bisect import bisect_left
import pysam
import bx.bbi.bigwig_file
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
import multinomial_cov


import warnings
warnings.filterwarnings('error')


class SignalTrack(PA.Track):
    """Class for getting V-plot signal"""
    def __init__(self, chrom, start, end):
        PA.Track.__init__(self, chrom, start, end, "signal")
    def calculateSignal(self, mat, vmat):
        offset=self.start-mat.start-vmat.w
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate signal")
        self.vals = signal.correlate(mat.get(vmat.i_lower,vmat.i_upper,
                                              mat.start + offset, mat.end - offset),
                                       vmat.mat,mode = 'valid')[0]

class NormSignalTrack(PA.Track):
    """Class for storing normalized signal track"""
    def __init__(self, chrom, start, end):
        PA.Track.__init__(self, chrom, start, end, "normalized signal")
    def calculateNormSignal(self, raw, bias):
        self.vals = raw.get(self.start, self.end) - bias.get(self.start,self.end)

class BiasTrack(PA.Track):
    """Class for getting Bias Signal Track-- Background model"""
    def __init__(self, chrom, start, end):
        PA.Track.__init__(self, chrom, start, end, "bias")
    def calculateBackgroundSignal(self, mat, vmat, nuc_cov):
        offset=self.start-mat.start-vmat.w
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate signal")
        self.vmat = vmat
        self.bias_mat = mat
        self.cov = PA.CoverageTrack(self.chrom, self.start, self.end)
        self.cov.calculateCoverage(self.bias_mat, vmat.i_lower,
                                   vmat.i_upper, vmat.w*2+1)
        self.nuc_cov = nuc_cov.vals
        self.vals = signal.correlate(self.bias_mat.get(vmat.i_lower,vmat.i_upper,
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
        bias_mat = bias_mat.get(vmat.i_lower,vmat.i_upper,position - vmat.w,position  + vmat.w + 1)
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
        var = multinomial_cov.makeCov(self.probs, flatv, self.reads)
        return np.sqrt(var)
    def analMean(self):
        return np.sum(self.prob_mat * self.vmat.mat * self.reads)



def norm(x, v, w, mean):
    """compute values of normal pdf with given mean and sd at values in x"""
    norm = (1.0/(np.sqrt(2*np.pi*v)) *
        np.exp(-(x - mean)**2/(2*v)))
    norm = norm * (w/max(norm))
    return norm

class Nucleosome:
    """Class for storing information about a single nucleosome"""
    def __init__(self, pos,nuctrack):
        self.chrom = nuctrack.chrom
        self.pos = pos
        self.nfr_cov = nuctrack.nfr_cov.get(pos = pos)
        self.nuc_cov = nuctrack.nuc_cov.get(pos = pos)
        self.nuc_signal = nuctrack.nuc_signal.get(pos = pos)
        self.norm_signal = nuctrack.norm_signal.get(pos = pos)
        self.smoothed = nuctrack.smoothed.get(pos= pos)
    def getLR(self,nuctrack):
        mat = nuctrack.mat.get(nuctrack.params.lower,nuctrack.params.upper, 
                                self.pos - nuctrack.params.vmat.w, self.pos + nuctrack.params.vmat.w +1)
        null_mat = nuctrack.bias_mat.get(nuctrack.params.lower,nuctrack.params.upper,      
                                self.pos - nuctrack.params.vmat.w, self.pos + nuctrack.params.vmat.w +1)
        bias_mat =nuctrack.bias_mat_prenorm.get(nuctrack.params.lower,nuctrack.params.upper,
                                self.pos - nuctrack.params.vmat.w, self.pos + nuctrack.params.vmat.w +1)
        nuc_model = nuctrack.params.vmat.mat * bias_mat
        nuc_model = nuc_model / np.sum(nuc_model)
        null_model = null_mat / np.sum(null_mat)
        nuc_lik = np.sum(np.log(nuc_model) * mat)
        null_lik = np.sum(np.log(null_model) * mat)
        self.lr = nuc_lik - null_lik
    def getZScore(self, nuctrack):
        s = SignalDistribution(self.pos, nuctrack.params.vmat, nuctrack.bias_mat,
                                self.nuc_cov)
        std = s.analStd()
        self.z = self.norm_signal / std
    def getOcc(self, nuctrack):
        self.occ = nuctrack.occ.get(pos = self.pos)
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
        index = self.pos - nuctrack.start
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
    def writeNucPos(self, handle):
        out = [self.chrom, self.pos, self.pos + 1, self.z, self.occ, self.lr,
               self.norm_signal, self.smoothed, self.nuc_signal, self.nuc_cov, self.nfr_cov,
                self.fuzz]
        handle.write("\t".join(map(str, out)) + "\n")

class NFR:
    """Class for storing information about a single NFR region"""
    def __init__(self, left, right, nuctrack):
        self.chrom = nuctrack.chrom
        self.leftpos = left
        self.rightpos = right
        self.nuc = np.max(nuctrack.nuc_signal.get(left,right))
        self.nfr = np.mean(nuctrack.nfr_cov.get(left,right))
        self.occ = np.mean(nuctrack.occ.get(left,right))
        self.max_occ = np.max(nuctrack.occ.get(left,right))
    def writeNFRPos(self, handle):
        out = [self.chrom, self.leftpos, self.rightpos, self.occ, self.nfr, self.nuc, self.max_occ]
        handle.write("\t".join(map(str, out)) + "\n")


class NucParameters:
    """Class for storing parameters related to nucleosome calling"""
    def __init__(self, vmat, insertsizes, sd = 25, nonredundant_sep = 120, redundant_sep = 25,
                 min_z = 3, min_lr = 0, downsample = None, min_reads = 1,
                 seed = 500, min_nfr_len = 1, max_nfr_len = 1000, max_nfr_occ = 0.2, min_nfr_ins = 0.05,
                out = None, bias = None, gdna = None, write_all = False, occ_track = None, 
                bam = None, occ_params = None):
        self.vmat = vmat
        self.lower = vmat.i_lower
        self.upper= vmat.i_upper
        self.window = vmat.mat.shape[1]
        self.insertsizes= insertsizes
        self.min_reads = min_reads
        self.min_z = min_z
        self.min_lr = min_lr
        self.smooth_sd = sd
        self.redundant_sep = redundant_sep
        self.nonredundant_sep = nonredundant_sep
        self.downsample = downsample
        self.seed = seed
        self.min_nfr_len = min_nfr_len
        self.max_nfr_len = max_nfr_len
        self.max_nfr_occ = max_nfr_occ
        self.min_nfr_ins = min_nfr_ins
        self.out = out
        self.bias = bias
        self.gdna = gdna
        self.write_all = write_all
        self.occ_track = occ_track
        self.bam = bam
        self.occ_params = occ_params

class NucChunk(PA.ChromChunk):
    """Class for storing and determining collection of nucleosome positions
    """
    def __init__(self, chrom, start, end):
        PA.ChromChunk.__init__(self, chrom, start, end)
    def initialize(self, parameters):
        self.params = parameters
    def extract_reads(self, bamfile):
        self.reads = PA.ReadList(self.chrom, self.start - self.params.window,
                                 self.end + self.params.window)
        self.reads.extract_reads(bamfile,i_upper=self.params.upper, downsample = self.params.downsample,
                                 seed = self.params.seed)
    def getReadMat(self):
        self.mat = PA.ReadMat2D(self.chrom, self.start - self.params.window,
                                 self.end + self.params.window, 0, self.params.upper)
        self.mat.makeReadMat(self.reads)
    def makeBiasMat(self,  biasfile = None, gdna = False):
        self.bias_mat = PA.BiasMat2D(self.chrom, self.start - self.params.window,
                                 self.end + self.params.window, 0, self.params.upper)
        if biasfile and gdna:
            bias_track = PA.InsertionBiasTrack(self.chrom, self.start - self.params.window - self.params.upper/2,
                                  self.end + self.params.window + self.params.upper/2 + 1, log = False)
            bias_track.read_track(biasfile)
            self.bias_mat.makeBiasMat(bias_track)
        elif biasfile:
            bias_track = PA.InsertionBiasTrack(self.chrom, self.start - self.params.window - self.params.upper/2,
                                  self.end + self.params.window + self.params.upper/2 + 1, log = True)
            bias_track.read_track(biasfile)
            self.bias_mat.makeBiasMat(bias_track)
        self.bias_mat_prenorm = PA.BiasMat2D(self.chrom, self.start - self.params.window,
                                 self.end + self.params.window, 0, self.params.upper)
        self.bias_mat_prenorm.mat = copy(self.bias_mat.mat)
        self.bias_mat.normByInsertDist(self.params.insertsizes)
    def getNucSignal(self):
        """Gets Nucleosome Signal Track"""
        self.nuc_cov = PA.CoverageTrack(self.chrom, self.start,
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
        self.nfr_cov = PA.CoverageTrack(self.chrom, self.start, self.end)
        self.nfr_cov.calculateCoverage(self.mat, 0, self.params.lower,
                                                        self.params.window)
    def smoothSignal(self):
        """Smooth thenormalized signal track"""
        window_len = 6 * self.params.smooth_sd + 1
        self.smoothed = PA.Track(self.chrom,self.start,self.end, "Smooth Signal")
        tmp = copy(self.norm_signal.vals)
        self.smoothed.assign_track(tmp)
        self.smoothed.vals[ self.smoothed.vals < 0] = 0
        self.smoothed.smooth_track(window_len, window = "gaussian",
                             sd = self.params.smooth_sd, mode = 'same',
                             norm = True)
    def getOcc(self, occ_bwh = None):
        """gets occupancy track-- either reads in from bw handle given, or makes new"""
        if occ_bwh is None and self.params.occ_params is None:
            raise Exception("Error with getOcc functoin for NucChunk.  Need either bw handle or OccCalcParams")
        elif occ_bwh is not None and self.params.occ_params is not None:
            raise Warning("both occupancy track and parameters given-- will only use track")
        if occ_bwh is not None:
            self.occ = PA.Track(self.chrom,self.start,self.end,"Occupancy")
            self.occ.read_track(occ_bwh)
        else:
            self.occ = Occ.OccupancyTrack(self.chrom, self.start, self.end)
            self.occ.calculateOccupancyMLE_smooth(self.mat, self.params.occ_params)
    def findAllNucs(self):
        """Find peaks in data"""
        self.nuc_collection = {}
        combined = self.norm_signal.vals + self.smoothed.vals
        #find peaks in normalized sigal
        cands1 = PA.call_peaks(combined, min_signal = 0,
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
        self.nonredundant = PA.reduce_peaks( self.sorted_nuc_keys, 
                                            map(lambda x: self.nuc_collection[x].z, self.sorted_nuc_keys),
                                                self.params.nonredundant_sep)
        self.redundant = np.setdiff1d(self.sorted_nuc_keys, self.nonredundant)
    def fit(self):
        x = np.linspace(0,self.length -1,self.length)
        fit = np.zeros(self.length)
        for nuc in self.sorted_nuc_keys:
            self.nuc_collection[nuc].getFuzz(self)
            fit += norm(x,self.nuc_collection[nuc].fuzz**2, self.nuc_collection[nuc].weight, self.nuc_collection[nuc].fit_pos)
        self.fitted = PA.Track(self.chrom, self.start, self.end,
                            "Fitted Nucleosome Signal")
        self.fitted.assign_track(fit)
    def findNFRs(self):
        """Find NFR regions based on nfr reads, nucleosomes"""
        self.nfrs = {}
        nucpos = [i for i in self.sorted_nuc_keys if self.nuc_collection[i].occ > self.params.max_nfr_occ]
        gaps =[nucpos[i] - nucpos[i-1] for i in range(1,len(nucpos))]
        for i in range(len(gaps)):
            if gaps[i]<self.params.min_nfr_len+147 or gaps[i]>self.params.max_nfr_len+147:
                continue
            leftnucbound = nucpos[i] + 73
            rightnucbound = nucpos[i+1] - 73
            if max(self.occ.vals[leftnucbound:rightnucbound]) < self.params.max_nfr_occ and np.mean(self.ins.vals[leftnucbound:rightnucbound]) > self.params.min_nfr_ins:
                self.nfrs[leftnucbound] = NFR(leftnucbound + self.start, rightnucbound + self.start + 1,self)
    def makeInsertionTrack(self):
        """make insertion track for chunk"""
        self.ins = PA.InsertionTrack(self.chrom, self.start,self.end)
        self.ins.calculateInsertions(self.reads)
    def process(self, params, bamfile, occfile, biasfile = None, gdna = False):
        """wrapper to carry out all methods needed to call nucleosomes and nfrs"""
        self.initialize(params)
        self.extract_reads(bamfile)
        self.getReadMat()
        if biasfile is not None:
            self.makeBiasMat(biasfile)
        elif gdna:
            self.makeBiasMat(biasfile, gdna = True)
        else:
            self.makeBiasMat()
        self.getNucSignal()
        self.getNFR()
        self.smoothSignal()
        self.getOcc(occ_bwh = occfile)
        self.findAllNucs()
        self.fit()
        self.makeInsertionTrack()
        self.findNFRs()
    def removeData(self):
        """remove data from chunk-- deletes all attributes"""
        names = self.__dict__.keys()
        for name in names:
            delattr(self,name)

class NucChunkSet(PA.ChunkSet):
    """Class for storing sets of NucChunks
    """
    def __init__(self, name, chunks):
        PA.ChunkSet.__init__(self, name, chunks)
    def process(self, params):
        """process chunkset-- for all chunks, perform all methods needed to 
        obtain and write nucleosome signal and positions
        """
        if params.bias is not None and params.gdna is not None:
            raise Exception("NucChunkSet.process shoud only have either a bias or gdna input, not both")
        bamfile = pysam.Samfile(params.bam, "rb")
        if params.occ_track:
            occfile = bx.bbi.bigwig_file.BigWigFile(open(params.occ_track,'rb'))
        else:
            occfile = None
        if params.bias is not None:
            biasfile = bx.bbi.bigwig_file.BigWigFile(open(params.bias,'rb'))
        elif params.gdna is not None:
            biasfile = bx.bbi.bigwig_file.BigWigFile(open(params.bias,'rb'))
        else:
            biasfile = None
        handles = {}
        if params.write_all:
            outputs = ['nucpos','nucpos.redundant','nfrpos','nucsig','smoothsig','ins','nuc_cov',
                       'nfr_cov','background','rawsig','fitted']
        else:
            outputs = ['nucpos','nucpos.redundant','nfrpos','nucsig','smoothsig','ins']
        if params.occ_track is None:
            outputs.append('occ')
        for i in outputs:
            handles[i] = open('tmp_nucleoatac/' + self.name + '.' + params.out + '.' +
                            i + '.tmp.bed','w')
        for chunk in self.chunks:
            if params.bias:
                chunk.process(params, bamfile, occfile, biasfile, gdna = False)
            elif params.gdna:
                chunk.process(params, bamfile, occfile, biasfile, gdna = True)
            else:
                chunk.process(params, bamfile, occfile)
            for key in chunk.nonredundant:
                chunk.nuc_collection[key].writeNucPos(handles['nucpos'])
            for key in chunk.redundant:
                chunk.nuc_collection[key].writeNucPos(handles['nucpos.redundant'])
            for key in chunk.nfrs.keys():
                chunk.nfrs[key].writeNFRPos(handles['nfrpos'])
            chunk.ins.write_track(handles['ins'])
            chunk.norm_signal.write_track(handles['nucsig'])
            chunk.smoothed.write_track(handles['smoothsig'])
            if params.write_all:
                chunk.nuc_signal.write_track(handles['rawsig'])
                chunk.nuc_cov.vals = chunk.nuc_cov.vals / params.window * 10.0
                chunk.nuc_cov.write_track(handles['nuc_cov'])
                chunk.bias.write_track(handles['background'])
                chunk.nfr_cov.vals = chunk.nfr_cov.vals / params.window * 10.0
                chunk.nfr_cov.write_track(handles['nfr_cov'])
                chunk.fitted.write_track(handles['fitted'])
            if params.occ_track is None:
                chunk.occ.write_track(handles['occ'])
            chunk.removeData()
        for key in handles.keys():
            handles[key].close()
        bamfile.close()

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
            out.extend([NucChunk(chrom, start, end)])
        infile.close()
        return out


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
            sets.append(NucChunkSet(str(j).zfill(z),bed_list[j:(j+1)]))
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
            sets.append(NucChunkSet(str(j).zfill(z),bed_list[a:b]))
    return sets
