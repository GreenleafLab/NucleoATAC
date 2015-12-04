"""
General tools for dealing with ATAC-Seq data using Python.

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

from pkg_resources import resource_filename
import os
import numpy as np
import scipy.signal as signal
import pyatac.seq as seq
from pyatac.tracks import Track
from pyatac.utils import smooth


def pwm_parse(name):
    out = resource_filename('pyatac.pwm', name + '.PWM.txt')
    if os.path.isfile(out):
        return out
    else:
        return name


class PWM:
    """Class for storing Tn5 preference pwm"""
    def __init__(self, mat, up, down, nucleotides):
        self.mat = mat
        self.up = up
        self.down = down
        self.nucleotides = nucleotides
    def save(self,filename):
        """write text output description of PWM object attributes"""
        out=open(filename,'w')
        out.write('#PWM Descriptor File\n')
        out.write('#Contains PWM and pertinent information\n')
        out.write('#up\n')
        out.write(str(self.up)+'\n')
        out.write('#down\n')
        out.write(str(self.down)+'\n')
        out.write('#nucleotides\n')
        out.write("\t".join(self.nucleotides)+'\n')
        out.write('#mat\n')
        for row in self.mat:
            out.write("\t".join(map(str,row))+'\n')
        out.close()
    @staticmethod
    def open(name):
        """Create PWM object from text descriptor file"""
        filename = pwm_parse(name)
        infile = open(filename,'r')
        state = ''
        mat = []
        for line in infile:
            if '#up' in line:
                state = 'up'
            elif '#down' in line:
                state = 'down'
            elif '#mat' in line:
                state = 'mat'
            elif '#nucleotides' in line:
                state = 'nucleotides'
            elif state == 'up':
                up = int(line.strip('\n'))
            elif state == 'down':
                down = int(line.strip('\n'))
            elif state == 'nucleotides':
                nucleotides = line.strip('\n').split()
            elif state == 'mat':
                mat.append(map(float,line.strip('\n').split('\t')))
        infile.close()
        try:
            new = PWM(np.array(mat), up, down, nucleotides)
        except NameError:
            raise Exception("PWM decriptor file appeas to be missing some\
needed components")
        return new




class InsertionBiasTrack(Track):
    """Class for getting and storing insertion bias"""
    def __init__(self, chrom, start, end, log = True):
        Track.__init__(self, chrom, start, end, name = "insertion bias", log = log)
    def computeBias(self, fasta, chromDict, pwm):
        """compute bias track based on sequence and pwm"""
        self.slop(chromDict, up = pwm.up, down = pwm.down)
        sequence = seq.get_sequence(self, fasta)
        seqmat = seq.seq_to_mat(sequence, pwm.nucleotides)
        self.vals = signal.correlate(seqmat,np.log(pwm.mat),mode='valid')[0]
        self.start += pwm.up
        self.end -= pwm.down
    def get(self, start = None, end = None, pos = None, log = None):
        """Obtain values for particular region"""
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



