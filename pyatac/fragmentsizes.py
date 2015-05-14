"""
Classes for working with fragment distribution

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

import numpy as np
import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()})
from pyatac.fragments import getAllFragmentSizes, getFragmentSizesFromChunkList




class FragmentSizes:
    """Class for storing fragment size distribution"""
    def __init__(self, lower, upper, atac = True, vals = None):
        self.lower = lower
        self.upper = upper
        self.vals = vals
        self.atac = atac
    def calculateSizes(self, bamfile, chunks = None):
        if chunks is None:
            sizes = getAllFragmentSizes(bamfile, self.lower, self.upper, atac = self.atac)
        else:
            sizes = getFragmentSizesFromChunkList(chunks, bamfile, self.lower, self.upper, atac = self.atac)
        self.vals = sizes / (np.sum(sizes) + (np.sum(sizes)==0))
    def get(self, lower = None, upper = None, size = None):
        if size:
            try:
                return self.vals[size - self.lower]
            except:
                raise Exception("Looks like size doesn't match FragmentSizes")
        else:
            if lower is None:
                lower = self.lower
            if upper is None:
                upper = self.upper
            y1 = lower - self.lower
            y2 = upper - self.lower
            try:
                return self.vals[y1:y2]
            except:
                raise Exception("Looks like dimensions from get probaby don't match FragmentSizes")
    def save(self, filename):
        """Save Fragment Distribution information"""
        f = open(filename,"w")
        f.write("#lower\n")
        f.write(str(self.lower)+"\n")
        f.write("#upper\n")
        f.write(str(self.upper)+"\n")
        f.write("#sizes\n")
        f.write("\t".join(map(str,self.get()))+"\n")
        f.close()
    @staticmethod
    def open(filename):
        """Create FragmentDistribution object from text descriptor file"""
        infile = open(filename,'r')
        state = ''
        for line in infile:
            if '#lower' in line:
                state = 'lower'
            elif '#upper' in line:
                state = 'upper'
            elif '#sizes' in line:
                state = 'sizes'
            elif '#' in line:
                state = 'other'
            elif state == 'lower':
                lower = int(line.strip('\n'))
            elif state == 'upper':
                upper = int(line.strip('\n'))
            elif state == 'sizes':
                fragmentsizes = np.array(map(float,line.rstrip("\n").split("\t")))
        try:
            new = FragmentSizes(lower, upper, vals = fragmentsizes)
        except NameError:
            raise Exception("FragmentDistribution decriptor file appeas to be missing some\
needed components")
        infile.close()
        return new



