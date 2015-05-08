#!/usr/bin/env python
"""
General tools for dealing with ATAC-Seq data using Python.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import subprocess
import numpy as np
from scipy import signal
import pysam
import warnings
from copy import copy

#Run shell command
def shell_command(cmd):
    """Conduct shell command."""
    output = subprocess.check_output(cmd, shell = True)
    return(output)


#Smoothing function
def smooth(sig, window_len, window='flat', sd = None, mode = 'valid',
           norm = True):
    """smoothes input signal using either flat or gaussian window

    options for window are "flat" and "gaussian"
    if window is "gaussian", sd should be provided.
    for guassian default sd is (window_len-1)/6
    norm means whether window should integrate to 1
    """
    if window not in ['flat','gaussian']:
        raise Exception("Incorrect window input for smooth. Options are flat, gaussian")
    if window_len%2 != 1:
        warnings.warn("Window length is even number.  Needs to be odd so adding 1.")
        window_len += 1
    if window=='gaussian' and sd is None:
        sd = (window_len-1)/6.0
    if window=="gaussian":
        w = signal.gaussian(window_len,sd)
    if window=="flat":
        w = np.ones(window_len)
    sig_nonan = copy(sig)
    sig_nonan[np.isnan(sig)] = 0
    smoothed = np.convolve(w, sig_nonan, mode = mode) 
    if norm:
        norm_sig = np.ones(len(sig))
        norm_sig[np.isnan(sig)] = 0
        smoothed_norm = np.convolve(w, norm_sig, mode = mode)
        smoothed_norm[smoothed_norm==0] = np.nan
        smoothed = smoothed / smoothed_norm
    return smoothed



def reduce_peaks(peaks,sig, sep):
    """Greedy algorithm for taking peaks and turning to set with at least sep distance
        between peaks.  First include peak with highest sig value, then next greatest
        not within sep distance from that one, and so on"""
    exclude = np.zeros(peaks.size)
    keep = np.zeros(peaks.size)
    st = np.argsort(sig)
    j=peaks.size-1
    while j >=0:
        ind = st[j]
        j+= -1
        if exclude[ind]==0:
            keep[ind] = 1
            exclude[ind]=1
            k = ind - 1
            while k>=0 and  (peaks[ind] - peaks[k]) < sep:
                exclude[k]=1
                k += -1
            k = ind + 1
            while k < peaks.size and (peaks[k]-peaks[ind]) < sep:
                exclude[k]=1
                k += 1
    return peaks[keep ==1]



def call_peaks(sigvals, min_signal = 0, sep = 120, boundary = None, order =1):
    """Greedy algorithm for peak calling-- first call all local maxima,
    then call greatest maxima as a peak, then next greatest that isn't within
    'sep' distance of that peak, and so on"""
    if sum(np.isnan(sigvals))>0:
        if sum(np.isnan(sigvals))==len(sigvals):
            return np.array([])
        else:
            replace = min(sigvals[~np.isnan(sigvals)])
            sigvals[np.isnan(sigvals)]=replace
    if boundary is None:
        boundary = sep/2
    random = np.random.RandomState(seed = 25)
    l = len(sigvals)
    peaks = signal.argrelmax(sigvals *(1+
                        random.uniform(0,10**-12,l)),order=order)[0]
    peaks = peaks[sigvals[peaks] >= min_signal ]
    peaks = peaks[ peaks >= boundary ]
    peaks = peaks[ peaks < (l - boundary)]
    sig = sigvals[peaks]
    return reduce_peaks(peaks, sig, sep)

def read_chrom_sizes_from_fasta(fastafile):
    """get chromosome size information from fasta file"""
    out = {}
    fasta = pysam.FastaFile(fastafile)
    chr_names = fasta.references
    chr_lengths = fasta.lengths
    fasta.close()
    for i in range(len(chr_lengths)):
        out[chr_names[i]]=int(chr_lengths[i])
    return out


def read_chrom_sizes_from_bam(bamfile):
    """get chromosome size information from bamfile"""
    out = {}
    bam = pysam.Samfile(bamfile, "rb")
    chr_lengths=bam.lengths
    chr_names=bam.references
    bam.close()
    for i in range(len(chr_lengths)):
        out[chr_names[i]]=int(chr_lengths[i])
    return out

def read_chrom_sizes(sizesFile):
    """get chromosome size information from chromosome sizes file"""
    out = {}
    f = open(sizesFile,'r')
    for line in f:
        keys = line.rstrip("\n").split("\t")
        out[keys[0]] = int(keys[1])
    return out







