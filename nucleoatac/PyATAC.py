#!/usr/bin/env python
"""
General tools for dealing with ATAC-Seq data using Python.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import os
from scipy import signal
import numpy as np
import warnings
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pysam

def shell_command(cmd):
    """Conduct shell command.  If error code, raise error"""
    code = os.system(cmd)
    if code!=0:
        raise Exception("error with shell command: "+cmd)



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
    if norm:
        w = w/sum(w)
    smoothed = np.convolve(w, sig, mode = mode)
    return smoothed


#Basic classes for storing ATAC-seq reads at different levels
class Read:
    """Class that stores read position and length"""
    def __init__(self, chrom, left, insert):
        self.chrom = chrom
        self.left = left#0-based
        self.right = left + insert#1-based
        self.insert = insert #Physical length between insertion (inclusive)

class ReadList:
    """Class that stores list of reads"""
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.reads = []
    def addRead(self, read):
        """Add one read to thecollection"""
        self.reads.extend([read])
    def extract_reads(self, bamfile,i_upper = 2000, downsample = None, seed =500, atac = True):
        """Get reads that overlap with chunk

        bamfile: pysam object
        flank: region on each side of chunk boundary to also include
        downsample: a float between 0 and 1, None for no downsampling
        """
        if downsample is not None and downsample!=1:
            random = np.random.RandomState(seed)
        for p2_rds in bamfile.fetch(self.chrom, max(0, self.start - i_upper),
                                    self.end + i_upper):
            if p2_rds.mapq>=30 and p2_rds.is_proper_pair and not p2_rds.is_reverse:
                if atac:
                    #get left position
                    l_pos = p2_rds.pos + 4
                    #get insert size
                    #correct by 8 base pairs to be inserion to insertion
                    if p2_rds.tlen<0 and p2_rds.tlen==-(p2_rds.rlen):
                        ilen = abs(p2_rds.tlen) - 8
                    elif p2_rds.tlen<0:
                        continue
                    else:
                        ilen=p2_rds.tlen-8
                else:
                    l_pos = p2_rds.pos
                    if p2_rds.tlen<0 and p2_rds.tlen==-(p2_rds.rlen):
                        ilen = abs(p2_rds.tlen)
                    elif p2_rds.tlen<0:
                        continue
                    else:
                        ilen=p2_rds.tlen
                if downsample is None or downsample == 1:
                    self.addRead(Read(self.chrom, l_pos, ilen))
                else:
                    x = random.uniform(0,1)
                    if x <= downsample:
                        self.addRead(Read(self.chrom, l_pos, ilen))

class ChunkMat2D:
    """Class that stores reads in 2D matrix according to
    center position and length"""
    def __init__(self, chrom, start, end, i_lower, i_upper, mode = "centers"):
        modes = ["centers","ends"]
        if mode not in modes:
            raise Exception("Invalid mode. Expected one of: %s" % modes)
        self.mode = mode
        self.chrom = chrom
        self.i_lower = i_lower
        self.i_upper = i_upper
        self.start = start
        self.end = end
        self.ncol = end - start
        self.nrow = i_upper - i_lower
        self.mat = np.zeros((self.nrow, self.ncol))
    def get(self, i_lower = None, i_upper = None, start = None, end = None, flip = False):
        """get a subset (or all) of the matrix stored by ChunkMat2D object based on
        chromosal position and/or insert size"""
        if i_lower is None:
            i_lower = self.i_lower
        if i_upper is None:
            i_upper = self.i_upper
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        y1 = i_lower - self.i_lower
        y2 = i_upper - self.i_lower
        x1 = start - self.start
        x2 = end - self.start
        if not flip:
            try:
                return self.mat[y1:y2,x1:x2]
            except:
                raise Exception("Looks like dimensions from get probaby don't match ReadMat")
        else:
            if x1 < 1 or x2 > self.mat.shape[1] or y1 < 0 or y2 > self.mat.shape[0]:
                raise Exception("Looks like dimensions from get probaby don't match ReadMat")
            ncol = x2-x1
            if ncol%2==0:
                raise Exception("Can only flip mat if the width is odd!")
            else:
                new = np.zeros((y2-y1,ncol))
                for j in range(y1,y2):
                    if (j+self.i_lower)%2==1:
                        new[j,:] = self.mat[j,x1:x2][::-1]
                    else:
                        new[j,:] = self.mat[j,(x1-1):(x2)][::-1][1:]
                return new
    def assign(self, mat):
        """assign a matrix to object"""
        if mat.shape != self.mat.shape:
            raise Exception("Dimensions of input mat are wrong.  Uh oh!")
        self.mat = mat
    def save(self, filename):
        """Save object in a text file"""
        head = ",".join(map(str,[self.chrom,self.start,self.end,self.i_lower,self.i_upper]))
        np.savetxt(filename,self.mat,delimiter="\t", header = head)
    @staticmethod
    def open(filename):
        f = open(filename,'r')
        header = f.readline()
        f.close()
        elements = header.rstrip('\n').lstrip("#").split(',')
        mat = np.loadtxt(filename, skiprows=1)
        new= ChunkMat2D(elements[0],elements[1],elements[2],elements[3])
        new.assign(mat)
        return new
    def getIns(self):
        """Collape matrix into insertions.  Will reduce span on chromosome
        if mode is centers"""
        if self.mode == "centers":
            pattern = np.zeros((self.i_upper-self.i_lower,self.i_upper + (self.i_upper-1)%2))
            mid = self.i_upper/2
            for i in range(self.i_lower,self.i_upper):
                pattern[i-self.i_lower,mid+(i-1)/2]=1
                pattern[i-self.i_lower,mid-(i/2)]=1
            ins = signal.correlate2d(self.mat,pattern,mode="valid")[0]
            insertion_track = InsertionTrack(self.chrom,self.start + pattern.shape[1]/2, self.end - (pattern.shape[1]/2))
            insertion_track.assign_track(ins)
            return insertion_track
        elif self.mode == "ends":
            insertion_track = InsertionTrack(self.chrom,self.start, self.end)
            insertion_track.assign_track(np.sum(self.mat),axis=0)
            return insertion_track
    def plot(self, filename = None, title = None, i_lower = None,
             i_upper = None):
        """Plot 2d ReadMat"""
        if i_upper is None:
            i_upper = self.i_upper
        if i_lower is None:
            i_lower = self.i_lower
        fig = plt.figure()
        plt.imshow(self.get(i_lower= i_lower, i_upper = i_upper),
                   origin="lower",interpolation='nearest',
                extent=[self.start,self.end-1,i_lower,i_upper-1],cmap=cm.get_cmap('Greys'))
        plt.xlabel(self.chrom)
        plt.ylabel("Insert size")
        if title:
            plt.title(title)
        #plt.colorbar(shrink=0.8)
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            #Also save text output!
            filename2 = ".".join(filename.split(".")[:-1]+['txt'])
            np.savetxt(filename2,self.mat,delimiter="\t")
        else:
            fig.show()




class ReadMat2D(ChunkMat2D):
    """Class that stores read information in 2D matrix according"""
    def __init__(self, chrom, start, end, i_lower, i_upper, mode = "centers"):
        ChunkMat2D.__init__(self,chrom, start, end, i_lower, i_upper ,mode)
    def updateMat(self, read):
        row = read.insert - self.i_lower
        if self.mode == "centers":
            col = (read.insert-1)/2 + read.left - self.start
            if col>=0 and col<self.ncol and row<self.nrow and row>=0:
                self.mat[row, col] += 1
        else:
            col1 = read.left
            col2 = read.right - 1
            if col1>=0 and col1<self.ncol and row<self.nrow and row>=0:
                self.mat[row, col1] += 1
            if col2>=0 and col2<self.ncol and row<self.nrow and row>=0:
                self.mat[row, col2] += 1
    def makeReadMat(self, readlist):
        """Make 2D matrix using reads list"""
        for read in readlist.reads:
                self.updateMat(read)


class BiasMat2D(ChunkMat2D):
    """Class that stores read information in 2D matrix according"""
    def __init__(self, chrom, start, end, i_lower, i_upper, mode = "centers"):
        ChunkMat2D.__init__(self,chrom, start, end, i_lower, i_upper ,mode)
        self.mat = np.ones(self.mat.shape)
    def makeBiasMat(self, bias_track):
        """Make 2D matrix representing sequence bias preferences"""
        offset = self.i_upper/2
        bias = bias_track.get(self.start-offset,self.end+offset)
        if not bias_track.log:
            nonzero = np.where(bias !=0)[0]
            bias = np.log(bias + min(bias[nonzero]))
        pattern = np.zeros((self.i_upper-self.i_lower,self.i_upper + (self.i_upper-1)%2))
        mid = self.i_upper/2
        for i in range(self.i_lower,self.i_upper):
            pattern[i-self.i_lower,mid+(i-1)/2]=1
            pattern[i-self.i_lower,mid-(i/2)]=1
        for i in range(self.i_upper-self.i_lower):
            self.mat[i]=np.exp(np.convolve(bias,pattern[i,:],mode='valid'))
    def normByInsertDist(self, insertsizes):
        inserts = insertsizes.get(self.i_lower,self.i_upper)
        self.mat = self.mat * np.reshape(np.tile(inserts,self.mat.shape[1]),self.mat.shape,order="F")



class ChromChunk:
    """Class that stores reads for a particular chunk of the genome"""
    def __init__(self, chrom, start, end, weight = None, strand = None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = end - start
        self.weight = weight
        self.strand= strand


class ChunkSet:
    """Class for storing group of chunks"""
    def __init__(self, name, chunks):
        self.name = name
        self.chunks = chunks



class Chromosomes:
    """Class for storing chromosome information"""
    def __init__(self):
        self.chrs = {}
    def get_from_bam(self, bamfile):
        """Get chromosome info from bam"""
        bam = pysam.Samfile(bamfile, "rb")
        chr_lengths=bam.lengths
        chr_names=bam.references
        bam.close()
        for i in range(len(chr_lengths)):
            if chr_names[i]!='chrM':
                self.chrs[chr_names[i]]=int(chr_lengths[i])
    def write_chrs(self, filename):
        """write chromosome information to text file"""
        genome=open(filename,'w')
        for c in self.chrs.keys():
            genome.write('\t'.join(map(str,[c,self.chrs[c]]))+'\n')
        genome.close()



##Funcion for reading in bed file
def read_bed(bedfile, weight_col=None, strand_col = None):
    """Make a dictionary of Chromosome insances from a bedfile"""
    infile = open(bedfile,"r")
    out = []
    weight = None
    strand = None
    while 1:
        in_line = infile.readline().rstrip('\n').split("\t")
        if not in_line[0]:
            break
        if weight_col:
            weight=in_line[weight_col-1]
        if strand_col:
            strand = in_line[strand_col-1]
        out.extend([ChromChunk(in_line[0],int(in_line[1]), int(in_line[2]),
                                   weight, strand)])
    infile.close()
    return out


def makeChunkSetList(bed_list, number):
    """turn list of ChromChunks into a list of ChunkSet (sets of ChromChunks)"""
    l = len(bed_list)
    x = min(number,l)
    n = l/x
    sets = []
    z = int(np.log10(x)+1)
    for j in range(x):
        a = j*n
        if j == x-1:
            b=l
        else:
            b = (j+1)*n
        sets.append(ChunkSet(str(j).zfill(z),bed_list[a:b]))
    return(sets)


#Family of classes for storing Signal Tracks

class Track:
    """Generic class for various types of signal tracks"""
    def __init__(self, chrom, start, end, name = "track", vals=None , log = False):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = end - start
        self.log = log
    def assign_track(self, vals, start = None, end = None):
        """Assign values to track"""
        if start:
            self.start = start
        if end:
            self.end = end
        if len(vals)!= self.end - self.start:
            raise Exception("The values being assigned to track do not \
                    span the start to end of the track")
        self.vals = vals
    def write_track(self, handle, start = None, end = None, vals = None):
        """Write track to output file handle

        If vals are specified use those values
        Othersise use self.vals
        """
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        if vals is None:
            vals=self.vals
        if len(vals)!=self.end-self.start:
            print len(vals),self.end-self.start
            raise Exception("Error! Inconsistency between length of \
            values and start/end values")
        for k in range(start, end):
            k2 = k - start
            handle.write('\t'.join(map(str, [self.chrom, k, k + 1,
                                             vals[k2]])) + '\n')
    def read_track(self, bwh, start = None, end = None, flank = None):
        """Read track values from BigWig file handle"""
        if start:
            self.start = start
        if end:
            self.end = end
        if flank:
            self.start = self.start - flank
            self.end = self.end + flank
        self.vals = bwh.get_as_array(self.chrom, self.start,
                                     self.end)
        if self.vals is None:
            raise Exception("No values were read for desired track. \
            Perhaps the chromsome is not included in the bigwig file?")
    def log(self, pseudo = 1):
        if self.log:
            print "Logging a track that is already log..."
        adjusted = self.vals + pseudo
        self.vals = np.log(adjusted)
        self.log = True
    def exp(self):
        if not self.log:
            print "taking exponent of a non-logged track..."
        self.vals = np.exp(self.vals)
        self.log = False
    def smooth_track(self,  window_len, window='flat', sd = None,
                     mode = 'valid', norm = True):
        """ smoothing of track"""
        self.smoothed=True
        self.vals = smooth(self.vals, window_len, window = window, sd = sd,
                           mode = mode, norm = norm)
        if mode == 'valid':
            self.start = self.start + window_len/2
            self.end = self.end - window_len/2
            self.length = self.end - self.start
    def get(self, start = None, end = None, pos = None):
        if pos:
            try:
                return self.vals[pos-self.start]
            except:
                raise Exception("Looks like position given doesn't match track")
        else:
            if start is None:
                start = self.start
            if end is None:
                end = self.end
            x1 = start - self.start
            x2 = end - self.start
            try:
                return self.vals[x1:x2]
            except:
                raise Exception("Looks like dimensions from get probaby don't match track, or there are no vals in track")
    def plot(self, name = None, start = None, end = None, filename=None):
        """plot the values """
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        if name is None:
            name = self.name
        values = self.get(start,end)
        fig = plt.figure()
        plt.plot(range(start,end),values)
        plt.xlabel(self.chrom)
        plt.ylabel(name)
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            #Also save text output!
            filename2 = ".".join(filename.split(".")[:-1]+['txt'])
            out = np.vstack((range(start,end),values))
            np.savetxt(filename2,out,delimiter="\t")
        else:
            fig.show()

class InsertionBiasTrack(Track):
    """Class for getting and storing insertion bias"""
    def __init__(self, chrom, start, end, log = True):
        Track.__init__(self, chrom, start, end, name = "insertion bias", log = log)
    def get(self, start = None, end = None, pos = None, log = None):
        out = Track.get(self, start, end, pos)
        if log is None:
            return out
        elif log:
            if self.log:
                return out
            else:
                return np.log(out)
        else:
            if self.log:
                return np.exp(out)
            else:
                return out
    def getLocalbias(self, windowlen = 55, window = "flat"):
        """Make bias a reflection of score in a given window"""
        if self.log:
            ebias = np.exp(self.vals)
        else:
            ebias = self.vals
        smoothed = smooth(ebias,windowlen, window, norm = False)
        flank = windowlen/2
        if self.log:
            self.vals = np.log(ebias[flank:-flank]/(smoothed-ebias[flank:-flank]))
        else:
            self.vals = ebias[flank:-flank]/(smoothed-ebias[flank:-flank])
        self.start = self.start + flank
        self.end = self.end - flank

class InsertionTrack(Track):
    """Class for getting and storing insertion positions"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "insertions")
    def calculateInsertions(self, readlist, flank = 0, i_lower = 0, i_upper =2000):
        self.start = self.start - flank
        self.end = self.end + flank
        self.vals = np.zeros(self.end - self.start)
        for read in readlist.reads:
            if read.insert >=i_lower and read.insert<i_upper:
                left = read.left - self.start
                right = read.right - 1 - self.start
                if left>=0 and left<self.end - self.start:
                    self.vals[left] +=1
                if right>=0 and right<self.end - self.start:
                    self.vals[right] +=1


class CoverageTrack(Track):
    """Class for computing read center converage"""
    def __init__(self, chrom, start, end):
        Track.__init__(self, chrom, start, end, "coverage")
    def calculateCoverage(self,mat,i_lower,i_upper,window_len):
        offset=self.start-mat.start-(window_len/2)
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate coverage with desired window")
        i_lower=i_lower-mat.i_lower
        i_upper=i_upper-mat.i_lower
        if offset!=0:
            collapsed = np.sum(mat.mat[i_lower:i_upper,offset:-offset],axis=0)
        else:
            collapsed = np.sum(mat.mat[i_lower:i_upper,],axis=0)
        self.vals = smooth(collapsed, window_len, window="flat",
                            mode='valid',norm=False)
    def calculateCoverageSmooth(self,mat,i_lower,i_upper,window_len,sd):
        offset=self.start-mat.start-(window_len/2)
        if offset<0:
            raise Exception("Insufficient flanking region on \
                    mat to calculate coverage with desired window")
        i_lower=i_lower-mat.i_lower
        i_upper=i_upper-mat.i_lower
        if offset!=0:
            collapsed = np.sum(mat.mat[i_lower:i_upper,offset:-offset],axis=0)
        else:
            collapsed = np.sum(mat.mat[i_lower:i_upper,],axis=0)
        self.vals = smooth(collapsed, window_len, sd= sd, window="gaussian",
                            mode='valid',norm=False)




class InsertSizes:
    """Class for storing insert sze distribution"""
    def __init__(self, dist, i_lower, i_upper):
        self.i_lower = i_lower
        self.i_upper = i_upper
        self.vals = dist
    def get(self, i_lower = None, i_upper = None, size = None):
        if size:
            try:
                return self.vals[size - self.i_lower]
            except:
                raise Exception("Looks like size doesn't match InsertSizes")
        else:
            if i_lower is None:
                i_lower = self.i_lower
            if i_upper is None:
                i_upper = self.i_upper
            y1 = i_lower - self.i_lower
            y2 = i_upper - self.i_lower
            try:
                return self.vals[y1:y2]
            except:
                raise Exception("Looks like dimensions from get probaby don't match InsertSizes")
    def save(self, filename):
        """Save Insert Distribution information"""
        f = open(filename,"w")
        f.write("#i_lower\n")
        f.write(str(self.i_lower)+"\n")
        f.write("#i_upper\n")
        f.write(str(self.i_upper)+"\n")
        f.write("#isizes\n")
        f.write("\t".join(map(str,self.get()))+"\n")
        f.close()
    @staticmethod
    def open(filename):
        """Create InsertDistribution object from text descriptor file"""
        infile = open(filename,'r')
        state = ''
        for line in infile:
            if '#i_lower' in line:
                state = 'i_lower'
            elif '#i_upper' in line:
                state = 'i_upper'
            elif '#isizes' in line:
                state = 'isizes'
            elif '#' in line:
                state = 'other'
            elif state == 'i_lower':
                i_lower = int(line.strip('\n'))
            elif state == 'i_upper':
                i_upper = int(line.strip('\n'))
            elif state == 'isizes':
                insertsizes = np.array(map(float,line.rstrip("\n").split("\t")))
        try:
            new = InsertSizes(insertsizes,i_lower,i_upper)
        except NameError:
            raise Exception("InsertDistribution decriptor file appeas to be missing some\
needed components")
        infile.close()
        return new

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


def combine_tmp_bed(basename, track, directory = 'tmp_nucleoatac'):
    """Combines temporary bed files"""
    tmp_names = directory + '/*.'+ basename + '.'+ track +'.tmp.bed'
    bed_name = basename + '.' + track + '.bed'
    shell_command('cat ' + tmp_names + ' > ' + bed_name)
    shell_command('rm ' + tmp_names)

def combine_tmp_bed_to_bw(basename, track, directory = 'tmp_nucleoatac'):
    """combines temporary bed file and turns into bigwig"""
    tmp_names = directory + '/*.'+ basename + '.'+ track +'.tmp.bed'
    bed_name = directory + '/' + basename + '.' + track + '.bed'
    bw_name = basename + '.' + track + '.bw'
    genome = directory + '/genome.txt'
    shell_command('cat ' + tmp_names + ' > ' + bed_name)
    shell_command('rm ' + tmp_names)
    shell_command('bedGraphToBigWig '+ bed_name + ' ' + genome + ' ' + bw_name)
    shell_command('rm ' + bed_name)



