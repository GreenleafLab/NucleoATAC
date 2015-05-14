import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pyatac.tracks import InsertionTrack
import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
from fragments import makeFragmentMat

class ChunkMat2D:
    """Class that stores fragments in 2D matrix according to
    center position and length"""
    def __init__(self, chrom, start, end, lower, upper):
        self.chrom = chrom
        self.lower = lower
        self.upper = upper
        self.start = start
        self.end = end
        self.ncol = end - start
        self.nrow = upper - lower
        self.mat = np.zeros((self.nrow, self.ncol))
    def get(self, lower = None, upper = None, start = None, end = None, flip = False):
        """get a subset (or all) of the matrix stored by ChunkMat2D object based on
        chromosal position and/or insert size"""
        if lower is None:
            lower = self.lower
        if upper is None:
            upper = self.upper
        if start is None:
            start = self.start
        if end is None:
            end = self.end
        y1 = lower - self.lower
        y2 = upper - self.lower
        x1 = start - self.start
        x2 = end - self.start
        if not flip:
            try:
                return self.mat[y1:y2,x1:x2]
            except:
                raise Exception("Looks like dimensions from get probaby don't match Mat")
        else:
            if x1 < 1 or x2 > self.mat.shape[1] or y1 < 0 or y2 > self.mat.shape[0]:
                raise Exception("Looks like dimensions from get probaby don't match Mat")
            ncol = x2-x1
            if ncol%2==0:
                raise Exception("Can only flip mat if the width is odd!")
            else:
                new = np.zeros((y2-y1,ncol))
                for j in range(y1,y2):
                    if (j+self.lower)%2==1:
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
        head = ",".join(map(str,[self.chrom,self.start,self.end,self.lower,self.upper]))
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
        """Collape matrix into insertions.  Will reduce span on chromosome"""
        pattern = np.zeros((self.upper-self.lower,self.upper + (self.upper-1)%2))
        mid = self.upper/2
        for i in range(self.lower,self.upper):
            pattern[i-self.lower,mid+(i-1)/2]=1
            pattern[i-self.lower,mid-(i/2)]=1
        ins = signal.correlate2d(self.mat,pattern,mode="valid")[0]
        insertion_track = InsertionTrack(self.chrom,self.start + pattern.shape[1]/2, self.end - (pattern.shape[1]/2))
        insertion_track.assign_track(ins)
        return insertion_track
    def plot(self, filename = None, title = None, lower = None,
             upper = None):
        """Plot 2d ReadMat"""
        if upper is None:
            upper = self.upper
        if lower is None:
            lower = self.lower
        fig = plt.figure()
        plt.imshow(self.get(lower= lower, upper = upper),
                   origin="lower",interpolation='nearest',
                extent=[self.start,self.end-1,lower,upper-1],cmap=cm.get_cmap('Greys'))
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



class FragmentMat2D(ChunkMat2D):
    """Class that stores fragment information in 2D matrix according"""
    def __init__(self, chrom, start, end, lower, upper, atac = True):
        ChunkMat2D.__init__(self,chrom, start, end, lower, upper)
        self.atac = atac
    def updateMat(self, fragment):
        row = fragment.insert - self.lower
        if self.mode == "centers":
            col = (fragment.insert-1)/2 + fragment.left - self.start
            if col>=0 and col<self.ncol and row<self.nrow and row>=0:
                self.mat[row, col] += 1
        else:
            col1 = fragment.left
            col2 = fragment.right - 1
            if col1>=0 and col1<self.ncol and row<self.nrow and row>=0:
                self.mat[row, col1] += 1
            if col2>=0 and col2<self.ncol and row<self.nrow and row>=0:
                self.mat[row, col2] += 1
    def makeFragmentMat(self, bamfile):
        """Make 2D matrix"""
        self.mat = makeFragmentMat(bamfile, self.chrom, self.start, self.end, self.lower, self.upper, self.atac)


class BiasMat2D(ChunkMat2D):
    """Class that stores fragment probabilities in 2D matrix according to bias model"""
    def __init__(self, chrom, start, end, lower, upper):
        ChunkMat2D.__init__(self,chrom, start, end, lower, upper)
        self.mat = np.ones(self.mat.shape)
    def makeBiasMat(self, bias_track):
        """Make 2D matrix representing sequence bias preferences"""
        offset = self.upper/2
        bias = bias_track.get(self.start-offset,self.end+offset)
        if not bias_track.log:
            nonzero = np.where(bias !=0)[0]
            bias = np.log(bias + min(bias[nonzero]))
        pattern = np.zeros((self.upper-self.lower,self.upper + (self.upper-1)%2))
        mid = self.upper/2
        for i in range(self.lower,self.upper):
            pattern[i-self.lower,mid+(i-1)/2]=1
            pattern[i-self.lower,mid-(i/2)]=1
        for i in range(self.upper-self.lower):
            self.mat[i]=np.exp(np.convolve(bias,pattern[i,:],mode='valid'))
    def normByInsertDist(self, insertsizes):
        inserts = insertsizes.get(self.lower,self.upper)
        self.mat = self.mat * np.reshape(np.tile(inserts,self.mat.shape[1]),self.mat.shape,order="F")




