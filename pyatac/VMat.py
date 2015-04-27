#!/usr/bin/env python
"""
VMat Class

@author: Alicia Schep, Greenleaf Lab, Stanford University
"""

#Import necessary python modules
from scipy import signal, ndimage
import numpy as np
from copy import copy
import matplotlib.pyplot as plt

class VMat_Error(Exception):
    """Class for errors in VMat function"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class VMat:
    """Class for storing and processing V-plot matrix"""
    def __init__(self, mat, lower, upper):
        """
        Assumes Vplot is centered!
        Inputs:
        mat = matrix (as numpy array)
        lower = lower bound of insert sizes represented by mat
        upper = upper bound of insert sizes represented by mat
        """
        if mat.shape[0]!=upper-lower:
            raise VMat_Error("mat shape is not consistent with insert limits")
        self.mat = mat
        self.upper = upper
        self.lower = lower
        self.w = mat.shape[1]/2
    def trim(self,lower,upper,w):
        """reduce the size of the vplot

        lower is new lower bound
        upper is new upper bound
        w is new flanking region around center
        """
        up = upper-self.lower
        dn = lower-self.lower
        left = self.w - w
        right = self.w + w + 1
        if up > self.mat.shape[0] or dn < 0 or left < 0 or right > self.mat.shape[1]:
            raise VMat_Error("Mat is smaller than desired trim")
        self.mat = self.mat[dn:up,left:right]
        self.lower = lower
        self.upper = upper
        self.w = w
    def symmetrize(self):
        """Force the V-plot to be symmetric"""
        for j in range(self.lower,self.upper):
            i=j-self.lower
            if j%2==1:
                lefthalf = (self.mat[i,:(self.w+1)]+self.mat[i,self.w:][::-1])*0.5
                self.mat[i,:] = np.hstack((lefthalf,lefthalf[:-1][::-1]))
            else:
                righthalf = (self.mat[i,(self.w):-1]+self.mat[i,:self.w][::-1])*0.5
                self.mat[i,:] = np.hstack((righthalf[::-1],righthalf,self.mat[i,-1]))
    def flip(self, mode = 'same'):
        """Flip V-plot"""
        if mode == 'same':
            new = np.zeros(self.mat.shape)
            for j in range(self.lower,self.upper):
                i = j - self.lower
                if j%2==1:
                    new[i,:] = self.mat[i,][::-1]
                else:
                    new[i,:-1] = self.mat[i,:-1][::-1]
                    #for -1 postion don't actually have values
                    new[i,-1] = np.mean([self.mat[i,-1],self.mat[i,1]])
            self.mat = new
        elif mode == 'valid':
            new = np.zeros((self.mat.shape[0],self.mat.shape[1]-2))
            for j in range(self.lower,self.upper):
                i = j - self.lower
                if j%2==1:
                    new[i,:] = self.mat[i,1:-1][::-1]
                else:
                    new[i,:] = self.mat[i,:-1][::-1][1:]
            self.mat = new
            self.w += -1
        else:
            raise Exception("Mode must be one of 'same' or 'valid'")
    def smooth(self, sd = 1):
        """smooth v-plot using gaussian kernel"""
        self.mat = ndimage.filters.gaussian_filter(self.mat,sd,
                                                          mode='constant')
    def smooth1d(self, sd = 1, axis = 1):
        """smooth v-plot along one axis only"""
        self.mat = ndimage.filters.gaussian_filter1d(self.mat,sd,axis,
                                                          mode='nearest')
    def norm(self):
        """normalize v matrix so that signal minus even background will be 1 divided by base pairs in window"""
        tmp1 = self.mat / np.sum(self.mat)
        tmp2 = np.ones(self.mat.shape) * (1.0 / self.mat.size)
        self.mat = self.mat / (np.sum(self.mat * tmp1)- np.sum(self.mat * tmp2))
        self.mat = (self.mat / self.mat.shape[1]) * 10.0
    def norm_y(self,dist):
        """normalize vplot so insertsize matches supplied distribution"""
        for i in range(self.mat.shape[0]):
            self.mat[i] = self.mat[i]  * (dist.get(size = i + self.lower)/ np.sum(self.mat[i]))
    def converto1d(self):
        """convert the 2d matrix to a 1d representation of insertions"""
        self.one_d = np.zeros(self.upper + self.upper%2 +2*self.w+1)
        center = self.upper/2 + self.w
        for j in range(self.mat.shape[0]):
            for i in range(self.mat.shape[1]):
                ilen=j+self.lower
                val = copy(self.mat[j,i])
                if ilen%2==0:
                    self.one_d[center-(self.w-i)-(ilen/2)]+= val
                    self.one_d[center-(self.w-i)+(ilen/2)]+= val
                else:
                    self.one_d[center-(self.w-i)-(ilen/2)]+= val * 0.5
                    self.one_d[center-(self.w-i)+(ilen/2)]+= val * 0.5
                    self.one_d[center-(self.w-i)-(ilen/2+1)]+= val * 0.5
                    self.one_d[center-(self.w-i)+(ilen/2+1)]+= val * 0.5
        self.one_d = self.one_d / sum(self.one_d)
    def plot(self, mat=None, title=None, filename=None):
        """Plot current main matrix or specified matrix (of same dimensions)"""
        if mat is None:
            mat=self.mat
        elif mat.shape!=(self.upper-self.lower,self.w*2+1):
            raise VMat_Error("dimensions of input mat should match \
                                dim of vmat")
        fig = plt.figure()
        plt.imshow(mat,origin="lower",interpolation='nearest',
                extent=[-self.w,self.w,self.lower,self.upper-1])
        plt.xlabel("Position relative to dyad")
        plt.ylabel("Insert size")
        if title:
            plt.title(title)
        plt.colorbar(shrink=0.8)
        if filename:
            fig.savefig(filename)
            plt.close(fig)
        else:
            fig.show()
    def plot_1d(self,filename=None):
        """plot the 1d insertion representation of the matrix"""
        fig = plt.figure()
        xlim = len(self.one_d)/2
        plt.plot(range(-xlim,xlim+1),self.one_d)
        plt.vlines(-73,0,max(self.one_d)*1.1,linestyles='dashed')
        plt.vlines(73,0,max(self.one_d)*1.1,linestyles='dashed')
        plt.xlabel("Position relative to dyad")
        plt.ylabel("Insertion Frequency")
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            #Also save text output!
            filename2 = ".".join(filename.split(".")[:-1]+['txt'])
            np.savetxt(filename2,self.one_d,delimiter="\t")
        else:
            fig.show()
    def plot_insertsize(self,filename=None):
        """plot the insert size disribution in the main matrix"""
        fig = plt.figure()
        ins = np.sum(self.mat,axis=1)
        ins = ins/sum(ins)
        plt.plot(range(self.lower,self.upper),ins)
        plt.xlabel("Insert Size")
        plt.ylabel("Frequency")
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            #Also save text output!
            filename2 = ".".join(filename.split(".")[:-1]+['txt'])
            np.savetxt(filename2,ins,delimiter="\t")
        else:
            fig.show()
    def save(self,filename):
        """write text output description of VMat object attributes"""
        out=open(filename,'w')
        out.write('#VMat Descriptor File\n')
        out.write('#Contains VMat and pertinent information\n')
        out.write('#lower\n')
        out.write(str(self.lower)+'\n')
        out.write('#upper\n')
        out.write(str(self.upper)+'\n')
        out.write('#mat\n')
        for row in self.mat:
            out.write("\t".join(map(str,row))+'\n')
        out.close()
    @staticmethod
    def open(filename):
        """Create VMat object from text descriptor file"""
        infile = open(filename,'r')
        state = ''
        mat = []
        for line in infile:
            if '#lower' in line:
                state = 'lower'
            elif '#upper' in line:
                state = 'upper'
            elif '#mat' in line:
                state = 'mat'
            elif '#' in line:
                state = 'other'
            elif state == 'lower':
                lower = int(line.strip('\n'))
            elif state == 'upper':
                upper = int(line.strip('\n'))
            elif state == 'mat':
                mat.append(map(float,line.strip('\n').split('\t')))
        try:
            new = VMat(np.array(mat), lower, upper)
        except NameError:
            raise VMat_Error("VMat decriptor file appeas to be missing some\
needed components")
        infile.close()
        return new
