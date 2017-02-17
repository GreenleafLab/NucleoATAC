"""
Takes as input a raw vplot and processes it for use in nucleosome calling
Use insert distribution at nuclesomal regions in peaks for normalization

@author: Alicia Schep
"""

##### IMPORT MODULES #####
# import necessary for python
#import matplotlib as mpl
#mpl.use('PS')
import pyatac.VMat as V
from pyatac.fragmentsizes import FragmentSizes

def run_vprocess(args):
    """process vplot

    """
    vmat=V.VMat.open(args.vplot)
    #Trim, Symmetrize
    vmat.trim(args.lower,args.upper,args.flank)
    vmat.symmetrize()
    #insert size norm
    if args.sizes is not None:
     #read in fragmentsizes
        nuc_dist = FragmentSizes.open(args.sizes)
        vmat.norm_y(nuc_dist)
    ##Smooth
    if args.smooth > 0:
        vmat.smooth(sd=args.smooth)
    #normalize
    vmat.norm()
    #Make extra plots if requeted
    if args.plot_extra:
        vmat.autoCorr()
        vmat.plot_auto(args.out+'.vplot.Autocorr.eps')
        vmat.converto1d()
        vmat.plot_1d(args.out+'.vplot.InsertionProfile.eps')
        vmat.plot_insertsize(args.out+'.vplot.InsertSizes.eps')
    #make plot and save
    vmat.save(args.out+".VMat")
    vmat.plot(filename = args.out+".VMat.eps")



