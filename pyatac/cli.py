import argparse
from pyatac import __version__

def pyatac_main(args):
    """The Main function for calling pyatac

    """
    call = args.call
    if call == "vplot":
        from pyatac.make_vplot import make_vplot
        print('---------Making VPlot-----------------------------------------------------------')
        make_vplot(args)
    elif call == "bias_vplot":
        from pyatac.make_bias_vplot import make_bias_vplot
        if args.bg is None and args.fasta is None:
            print("\nMust input either --bg or --fasta! Exiting...")
            sys.exit(1)
        if args.bg is not None and args.fasta is not None:
            print("Warning-- both --bg and --fasta provided-- will us --bg for bias track and inrore --fasta")
        print('---------Making Bias VPlot------------------------------------------------------')
        make_bias_vplot(args)
    elif call == "signal":
        from pyatac.signal_around_sites import get_signal
        print('---------Getting signal around sites--------------------------------------------')
        get_signal(args)
    elif call == "ins":
        from pyatac.get_ins import get_ins
        print('---------Getting insertions to make track---------------------------------------')
        get_ins(args)
    elif call == "cov":
        from pyatac.get_cov import get_cov
        print('---------Getting insertions to make track---------------------------------------')
        get_cov(args)
    elif call == "nucleotide":
        from pyatac.get_nucleotide import get_nucleotide
        print ('---------Getting nucleotide content---------------------------------------')
        get_nucleotide(args)
    elif call == "bias":
        from pyatac.make_bias_track import make_bias_track
        print ('---------Making Tn5 Bias Track---------------------------------------')
        make_bias_track(args)
    elif call == "pwm":
        from pyatac.get_pwm import get_pwm
        print ('---------Making PWM from bam---------------------------------------')
        get_pwm(args)
    elif call == "sizes":
        from pyatac.get_sizes import get_sizes
        print ('---------Getting fragment sizes---------------------------------------')
        get_sizes(args)
    elif call == "counts":
        from pyatac.get_counts import get_counts
        print ('---------Getting fragment counts---------------------------------------')
        get_counts(args)





def pyatac_parser():
    """Prepares argparse object

    """
    argparser = argparse.ArgumentParser(description = "%(prog)s -- Utilitis for working with ATAC-Seq data",
                                        epilog = "For command line options for each command, type %(prog)s COMMAND -h")
    argparser.add_argument("--version", action="version", version="%(prog)s "+ __version__)
    subparsers = argparser.add_subparsers(dest = 'call' )

    add_nucleotide_parser( subparsers)

    add_vplot_parser( subparsers)

    add_bias_vplot_parser( subparsers)

    add_pwm_parser( subparsers)

    add_signal_parser( subparsers)

    add_ins_parser( subparsers)

    add_cov_parser( subparsers)

    add_bias_parser( subparsers)

    add_sizes_parser( subparsers)

    add_counts_parser( subparsers)

    return argparser

def add_counts_parser( subparsers):
    """Add argument parsers for the sizes utility

    """
    parser = subparsers.add_parser("counts", help = "pyatac function-- compute fragment counts within windows")
    group0 = parser.add_argument_group('Required', 'Necessary arguments')
    group0.add_argument('--bam', metavar='bam_file', help = 'Aligned reads', required=True)
    group0.add_argument('--bed', metavar='bed_file',required=True,
                        help = 'Windows in which to compute counts')
    group3 = parser.add_argument_group('Options', 'Optional settings')
    group3.add_argument('--out', metavar='output_basename', help = "Basename for output")
    group3.add_argument('--not_atac', action="store_false", dest = "atac", default = True, help = "Don't use atac offsets")
    group4 = parser.add_argument_group('Fragment size bounds option', 'Upper and lower limits')
    group4.add_argument('--lower',metavar="int", help = 'lower limit on insert \
        size. Default is 0',default=0,type=int)
    group4.add_argument('--upper', metavar="int",help="upper limit on insert \
        size.  Default is 500",default=500,type=int)
    return



def add_sizes_parser( subparsers):
    """Add argument parsers for the sizes utility

    """
    parser = subparsers.add_parser("sizes", help = "pyatac function-- compute fragment size distribution")
    group0 = parser.add_argument_group('Required', 'Necessary arguments')
    group0.add_argument('--bam', metavar='bam_file', help = 'Aligned reads', required=True)
    group2 = parser.add_argument_group('Find only fragmentsizes for regions of genome')
    group2.add_argument('--bed', metavar='bed_file',
                        help = 'Only compute size distribution for fragment centered within regions in bed file')
    group3 = parser.add_argument_group('Options', 'Optional settings')
    group3.add_argument('--out', metavar='output_basename', help = "Basename for output")
    group3.add_argument('--not_atac', action="store_false", dest = "atac", default = True, help = "Don't use atac offsets")
    group4 = parser.add_argument_group('Fragment sizs bounds option', 'Upper and lower limits')
    group4.add_argument('--lower',metavar="int", help = 'lower limit on insert \
        size. Default is 0',default=0,type=int)
    group4.add_argument('--upper', metavar="int",help="upper limit on insert \
        size.  Default is 500",default=500,type=int)
    group5 = parser.add_argument_group('Plot options', 'Make plots?')
    group5.add_argument('--no_plot', action="store_true", default=False,help = "Don't plot output")
    return

def add_bias_parser( subparsers):
    """Add argument parsers for the bias utility

    """
    parser = subparsers.add_parser("bias", help = "pyatac function-- compute Tn5 bias score")
    group0 = parser.add_argument_group('Required', 'Necessary arguments')
    group0.add_argument('--fasta', metavar = 'fasta_file', help = 'Accepts fasta file', required = True)
    group1 = parser.add_argument_group('PWM option', 'Designate PWM file or default will be used')
    group1.add_argument('--pwm', metavar = 'Tn5_PWM', help = "PWM descriptor file. Default is Human.PWM.txt included in package",
                         default = "Human")
    group2 = parser.add_argument_group('Find only bias for regions of genome')
    group2.add_argument('--bed', metavar='bed_file', help = 'Positions around which \
        to get nucleotide frequencies')
    group3 = parser.add_argument_group('Options', 'Optional settings')
    group3.add_argument('--out', metavar='output_basename', help = "Basename for output")
    group3.add_argument('--cores', metavar='int', help = "number of cores to use",
                        default = 1, type = int)
    return


def add_pwm_parser( subparsers):
    """Add argument parsers for the nucleotide utility

    """
    parser = subparsers.add_parser("pwm", help = "pyatac function-- get nucleotide (or di-nucleotide) content around sites")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--fasta', metavar = 'fasta_file', help = 'Accepts fasta file', required = True)
    group1.add_argument('--bam', metavar='bam_file', help = 'Reads around which to get nucleotide freq', required = True)
    group2 = parser.add_argument_group('Input sites','If bed provided, find nucleotide frequencies around read starts within bed region')
    group2.add_argument('--bed', metavar='bed_file', help = 'Regions from which to use reads')
    group3 = parser.add_argument_group('Options', 'Optional settings')
    group3.add_argument('--dinucleotide', action = 'store_true', help = "Compute dinucleotide frequencies instead of single nucleotide", default = False)
    group3.add_argument('--flank', metavar='int', help = "Bases away from insertion site to get frequencies for. Default is 10", default = 10, type = int)
    group3.add_argument('--lower',metavar="int", help = 'lower limit on insert size. default is 0', default=0,type=int)
    group3.add_argument('--upper', metavar="int", help="upper limit on insert size. default is 2000", default=2000,type=int)
    group3.add_argument('--not_atac', action="store_false", dest = "atac", default = True, help = "Don't use atac offsets")
    group3.add_argument('--no_sym', action="store_false", dest= "sym", default = True, help = "Don't symmetrize PWM")
    group3.add_argument('--out', metavar='output_basename', help = "Basename for output")
    group3.add_argument('--cores', metavar='int', help = "number of cores to use",
                        default = 1, type = int)
    return

def add_nucleotide_parser( subparsers):
    """Add argument parsers for the nucleotide utility

    """
    parser = subparsers.add_parser("nucleotide", help = "pyatac function-- get nucleotide (or di-nucleotide) content around sites")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--fasta', metavar = 'fasta_file', help = 'Accepts fasta file', required = True)
    group1.add_argument('--bed', metavar='bed_file', help = 'Positions around which \
        to get nucleotide frequencies', required = True)
    group3 = parser.add_argument_group('Options', 'Optional settings')
    group3.add_argument('--dinucleotide', action = 'store_true', help = "Compute dinucleotide frequencies instead of single nucleotide", default = False)
    group3.add_argument('--up', metavar='int', help = "Bases upstream of site to get frequencies for", default = 250, type = int)
    group3.add_argument('--down', metavar='int', help = "Bases downstream of site to get frequencies for", default = 250, type = int)
    group3.add_argument('--strand', metavar='int',help = "Column in bedfile with strand info (1-based)", type = int)
    group3.add_argument('--out', metavar='output_basename', help = "Basename for output")
    group3.add_argument('--cores', metavar='int', help = "number of cores to use",
                        default = 1, type = int)
    group3.add_argument('--norm', action = 'store_true', default = False, help = "Normalize by background frequencies")
    return




def add_vplot_parser( subparsers):
    """Add argument parsers for the vplot utility

    """
    parser = subparsers.add_parser("vplot", help = "pyatac function-- make vplot")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bed', metavar='bed_file' , help = 'Positions around which \
        to generate VPlot', required = True)
    group1.add_argument('--bam', metavar='bam_file', help = 'Accepts sorted BAM file',required=True)

    group2 = parser.add_argument_group('General options', '')
    group2.add_argument('--out', metavar='basename')
    group2.add_argument('--cores', metavar = 'int',default=1,
                    help='Number of cores to use',type=int)

    group3 = parser.add_argument_group('VMat option', 'Size, scaling of VPlot')
    group3.add_argument('--lower',metavar="int", help = 'lower limit on insert \
        size',default=0,type=int)
    group3.add_argument('--upper', metavar="int",help="upper limit on insert \
        size",default=250,type=int)
    group3.add_argument('--flank', metavar="int",
        help="how many bases on each side of site (or center of site) to include",
        type=int,default=250)
    group3.add_argument('--scale', action="store_true", default = False, help = "Scale each site")
    group3.add_argument('--weight',metavar="int",
        type=int,help="column in which weight information is included")
    group3.add_argument('--strand',metavar="int",
        type=int,help="column in which strand information is included")
    group3.add_argument('--not_atac', action="store_false", dest = "atac", default = True, help = "Don't use atac offsets")
    group5 = parser.add_argument_group('Plot options', 'Make plots?')
    group5.add_argument('--no_plot', action="store_true", default=False,help = "Don't plot output")
    group5.add_argument('--plot_extra', action="store_true",default=False,help = "Make some extra plots")
    return

def add_bias_vplot_parser( subparsers):
    """Add argument parsers for the bias vplot utility

    """
    parser = subparsers.add_parser("bias_vplot", help = "pyatac function-- make vplot")
    group0 = parser.add_argument_group('Required', 'Necessary arguments')
    group0.add_argument('--bed', metavar='bed_file' , help = 'Positions around which \
        to generate VPlot', required = True)
    group0.add_argument('--sizes', metavar='sizes_file', help = 'Accepts sizes file from pyatac sizes command',required=True)
    group1 = parser.add_argument_group('Bias Options',
                                       'Must either submit a tabix-indexed bedgraph file or a fasta file.')
    group1.add_argument('--bg', metavar='bias_file', help = 'Accepts tabix indexed file')
    group1.add_argument('--fasta', metavar='fasta_file', help = 'Accepts indexed fasta file')
    group1.add_argument('--pwm', metavar = 'Tn5_PWM', help = "PWM descriptor file. Default is Human.PWM.txt included in package",
                         default = "Human")
    group2 = parser.add_argument_group('General options', '')
    group2.add_argument('--out', metavar='basename')
    group2.add_argument('--cores', metavar = 'int',default=1,
                    help='Number of cores to use',type=int)
    group3 = parser.add_argument_group('VMat option', 'Size, scaling of VPlot')
    group3.add_argument('--lower',metavar="int", help = 'lower limit on insert \
        size',default=0,type=int)
    group3.add_argument('--upper', metavar="int",help="upper limit on insert \
        size",default=250,type=int)
    group3.add_argument('--flank', metavar="int",
        help="how many bases on each side of site (or center of site) to include",
        type=int,default=250)
    group3.add_argument('--scale', action="store_true", default = False)
    group3.add_argument('--weight',metavar="int",
        type=int,help="column in which weight information is included")
    group3.add_argument('--strand',metavar="int",
        type=int,help="column in which strand information is included")
    group4 = parser.add_argument_group('Plot options', 'Make plots?')
    group4.add_argument('--no_plot', action="store_true", default=False,help = "Don't plot output")
    group4.add_argument('--plot_extra', action="store_true",default=False,help = "Make some extra plots")
    return


def add_signal_parser( subparsers):
    """Add argument parsers for the signal utility

    """
    parser = subparsers.add_parser("signal", help = "pyatac function-- get signal around sites")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bed', metavar='bed_file' , help = 'Positions around which \
        to generate VPlot', required=True)
    group1.add_argument('--bg', metavar='bg_file', help = 'Accepts bedgraph file that is tabix indexed',required=True)
    group1.add_argument('--sizes', metavar ='genome_sizes_file', required = True, help = 'File with chromosome names in 1st col, sizes in 2nd')
    group2 = parser.add_argument_group('General options', '')
    group2.add_argument('--out', metavar='basename', help = 'basename for output')
    group2.add_argument('--cores', metavar = 'int', type = int,default=1,
                    help='Number of cores to use')
    group2.add_argument('--all', action='store_true', default=False,
                    help = "output csv file (gzipped) with signal track around\
                    all sites")
    group2.add_argument('--no_agg', action='store_true',default=False,
                    help = "Don't make a plot of aggregate or write up of aggregate")
    group3 = parser.add_argument_group('Bed options', 'Options related to bed intervals')
    group3.add_argument('--up', metavar="int",help="bases upstream of site to look",
        type=int,default=250)
    group3.add_argument('--down', metavar="int",help="bases dowstream site to look",
        type=int,default=250)
    group3.add_argument('--weight', metavar="int", type=int, default = None,
        help="Column with weight information. Signal for interval  will be weighted by value in column")
    group3.add_argument('--strand', metavar="int" ,type=int, default=None,
                    help = "Column in which strand information is included if strand is to be used")
    group4 = parser.add_argument_group('Signal options', 'Options related to signal')
    group4.add_argument('--exp', action='store_true', default=False,
                    help= "take exponent of value")
    group4.add_argument('--positive',action='store_true',default=False,
                   help = "Only include positive signal")
    group4.add_argument('--scale', action='store_true', default=False,
                    help = "scale each individual track by total\
                    signal value")
    group4.add_argument('--norm', action='store_true', default=False,
                    help = "normalize aggregate track by number of intervals")
    return

def add_ins_parser( subparsers):
    """Add argument parsers for the ins utility

    """
    parser = subparsers.add_parser("ins", help = "pyatac function-- get insertions")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bam', metavar='bam_file', help = 'Accepts sorted BAM file',required=True)
    group2 = parser.add_argument_group('General options', '')
    group2.add_argument('--bed', metavar='bed_file' , help = 'Regions in which to get insertions')
    group2.add_argument('--out', metavar='basename')
    group2.add_argument('--cores', metavar = 'int',default=1,
                    help='Number of cores to use',type=int)

    group3 = parser.add_argument_group('insertion option', 'Size range, smoothing')
    group3.add_argument('--lower',metavar="int", help = 'lower limit on insert \
        size',default=0,type=int)
    group3.add_argument('--upper', metavar="int",help="upper limit on insert \
        size",default=2000,type=int)
    group3.add_argument('--smooth', metavar = "int", type = int,
                         help = "smoothing window for guassian smoothing.  default is no smoothing")
    group3.add_argument('--not_atac', action="store_false", dest = "atac", default = True, help = "Don't use atac offsets")
    return

def add_cov_parser( subparsers):
    """Add argument parsers for the cov utility

    """
    parser = subparsers.add_parser("cov", help = "pyatac function-- get coverage")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bam', metavar='bam_file', help = 'Accepts sorted BAM file',required=True)
    group2 = parser.add_argument_group('General options', '')
    group2.add_argument('--bed', metavar='bed_file' , help = 'Regions in which to get insertions')
    group2.add_argument('--out', metavar='basename')
    group2.add_argument('--cores', metavar = 'int',default=1,
                    help='Number of cores to use',type=int)

    group3 = parser.add_argument_group('insertion option', 'Size range, smoothing')
    group3.add_argument('--lower',metavar="int", help = 'lower limit on insert \
        size',default=0,type=int)
    group3.add_argument('--upper', metavar="int",help="upper limit on insert \
        size",default=2000,type=int)
    group3.add_argument('--window', metavar = "int", type = int, default = 121,
                         help = "window for flat smoothing of coverage.  default is 121, should be odd")
    group3.add_argument('--scale', metavar = "float", type = float, default = 10,
                         help = "scaling value.  default is 10, corresponding to signal corresponding to # of fragment centers  per 10 bp. Use 1 for fragments per 1 bp.")
    group3.add_argument('--not_atac', action="store_false", dest = "atac", default = True, help = "Don't use atac offsets")
    return




