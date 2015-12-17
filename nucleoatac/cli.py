###-----------Import modules---------------####

import argparse
import nucleoatac.Magic
from nucleoatac import __version__

def nucleoatac_main(args):
    """The Main function for calling nucleoatac

    """
    #Parse options...
    call = args.call
    parser = nucleoatac_parser()
    if call == "occ":
        from nucleoatac.run_occ import run_occ
        print('---------Computing Occupancy and Nucleosomal Insert Distribution----------------')
        run_occ(args)
    elif call == "vprocess":
        from nucleoatac.run_vprocess import run_vprocess
        print('---------Processing VPlot-------------------------------------------------------')
        run_vprocess(args)
    elif call == "nuc":
        from nucleoatac.run_nuc import run_nuc
        print('---------Obtaining nucleosome signal and calling positions----------------------')
        run_nuc(args)
    elif call == "merge":
        from nucleoatac.merge import run_merge
        print('---------Merging----------------------------------------------------------------')
        run_merge(args)
    elif call == "nfr":
        from nucleoatac.run_nfr import run_nfr
        print('---------Calling NFR positions--------------------------------------------------')
        run_nfr(args)
    elif call == "run":
        occ_args = parser.parse_args(map(str,['occ','--bed',args.bed,'--bam',args.bam,
                                            '--fasta', args.fasta, '--pwm', args.pwm,
                                            '--out',args.out,'--cores',args.cores]))
        vprocess_args = parser.parse_args(['vprocess','--sizes',args.out+'.nuc_dist.txt','--out',args.out])
        nuc_args_list = ['nuc','--bed',args.bed,'--bam',args.bam,'--out',args.out,'--cores', str(args.cores),
                                        '--occ_track', args.out + '.occ.bedgraph.gz','--vmat', args.out + '.VMat',
                                        '--fasta', args.fasta, '--pwm', args.pwm, '--sizes', args.out + '.fragmentsizes.txt']
        if args.write_all:
            nuc_args_list.extend(['--write_all'])
        nuc_args = parser.parse_args(nuc_args_list)
        merge_args = parser.parse_args(['merge','--occpeaks',args.out +'.occpeaks.bed.gz','--nucpos',args.out+'.nucpos.bed.gz',
                                        '--out',args.out])
        nfr_args = parser.parse_args(['nfr','--bed', args.bed, '--occ_track', args.out + '.occ.bedgraph.gz', '--calls', 
                                        args.out + '.nucmap_combined.bed.gz','--out',args.out, '--fasta', args.fasta, 
                                        '--pwm', args.pwm , '--bam', args.bam])
        from nucleoatac.run_occ import run_occ
        from nucleoatac.run_vprocess import run_vprocess
        from nucleoatac.run_nuc import run_nuc
        from nucleoatac.merge import run_merge
        from nucleoatac.run_nfr import run_nfr
        print('---------Step1: Computing Occupancy and Nucleosomal Insert Distribution---------')
        run_occ(occ_args)
        print('---------Step2: Processing Vplot------------------------------------------------')
        run_vprocess(vprocess_args)
        print('---------Step3: Obtaining nucleosome signal and calling positions---------------')
        run_nuc(nuc_args)
        print('---------Step4: Making combined nucleosome position map ------------------------')
        run_merge(merge_args)
        print('---------Step5: Calling NFR positions-------------------------------------------')
        run_nfr(nfr_args)




def nucleoatac_parser():
    """Prepares argparse object

    """
    argparser = argparse.ArgumentParser(description = "%(prog)s -- Nucleosome Calling for ATAC-Seq",
                                        epilog = "For command line options for each command, type %(prog)s COMMAND -h")
    argparser.add_argument("--version", action="version", version="%(prog)s "+__version__)
    subparsers = argparser.add_subparsers(dest = 'call' )

    add_run_parser( subparsers)

    add_occ_parser( subparsers)

    add_vprocess_parser( subparsers)

    add_nuc_parser( subparsers)
    
    add_merge_parser( subparsers)

    add_nfr_parser( subparsers)

    return argparser

def add_occ_parser( subparsers):
    """Add argument parsers for the occ utility

    """
    parser = subparsers.add_parser("occ", help = "nucleoatac function:  Call nucleosome occupancy")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bed', metavar='bed_file' , help = 'Peaks in bed format', required=True)
    group1.add_argument('--bam', metavar='bam_file',
                    help = 'Sorted (and indexed) BAM file', required=True)
    group1.add_argument('--out', metavar='basename',
                    help="give output basename", required = True)
    group4 = parser.add_argument_group("Bias calculation information","Highly recommended. If fasta is not provided, will not calculate bias")
    group4.add_argument('--fasta', metavar = 'genome_seq',
                    help = 'Indexed fasta file')
    group4.add_argument('--pwm', metavar = 'Tn5_PWM', help = "PWM descriptor file. Default is Human.PWM.txt included in package", default = "Human")
    group2 = parser.add_argument_group('General Options', '')
    group2.add_argument('--sizes', metavar = 'fragmentsizes_file',
                       help = "File with fragment size distribution.  Use if don't want calculation of fragment size")
    group2.add_argument('--cores', metavar = 'int',default=1,
                    help='Number of cores to use',type=int)
    group3 = parser.add_argument_group('Occupancy parameter', 'Change with caution')
    group3.add_argument('--upper',metavar="int",default=251,
    help="upper limit in insert size. default is 251",type=int)
    group3.add_argument('--flank',metavar="int",default=60,
    help="Distance on each side of dyad to include for local occ calculation. Default is 60.",type=int)
    group3.add_argument('--min_occ', metavar = "float", default=0.1,type=float,
                    help="Occupancy cutoff for determining nucleosome distribution. Default is 0.1")
    group3.add_argument('--nuc_sep',metavar='int',default=120,type=int,
        help = "minimum separation between occupany peaks. Default is 120.")
    group3.add_argument('--confidence_interval', metavar='float',default=0.9,type=float,
        help = "confidence interval level for lower and upper bounds.  default is 0.9, should be between 0 and 1")
    group3.add_argument('--step', metavar = 'int', default = 5, type=int,
        help= "step size along genome for comuting occ. Default is 5.  Should be odd, or will be subtracted by 1")
    return

def add_merge_parser( subparsers):
    """Add argument parser for merge utility"""
    parser = subparsers.add_parser("merge", help = "nucleoatac function: Merge occ and nuc calls")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--occpeaks', metavar = 'occpeaks_file', required = True,
                        help = "Output from occ utility")
    group1.add_argument('--nucpos', metavar = 'nucpos_file', required = True,
                        help = "Output from nuc utility")
    group2 = parser.add_argument_group("Options","optional")
    group2.add_argument('--out', metavar = 'out_basename', help = "output file basename")
    group2.add_argument('--sep', metavar = 'min_separation', default = 120, 
                        help = "minimum separation between call")
    group2.add_argument('--min_occ', metavar = 'min_occ', default = 0.1,
                        help = "minimum lower bound occupancy of nucleosomes to be considered for excluding NFR. default is 0.1")

def add_nfr_parser( subparsers):
    """Add argument parser for nfr utility

    """
    parser = subparsers.add_parser("nfr", help = "nucleoatac function: Call NFRs")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bed', metavar='bed_file' , help = 'Peaks in bed format', required=True)
    group1.add_argument('--occ_track', metavar = 'occ_file', required = True,
                       help = "bgzip compressed, tabix-indexed bedgraph file with occcupancy track.")
    group1.add_argument('--calls', metavar = 'nucpos_file', required = True,
                       help = "bed file with nucleosome center calls")
    group6 = parser.add_argument_group("Insertion track options","Either input insertion track or bamfile")
    group6.add_argument('--ins_track', metavar = 'ins_file', 
                        help = "bgzip compressed, tabix-indexed bedgraph file with insertion track. will be generated if not included")
    group6.add_argument('--bam', metavar='bam_file',
                    help = 'Sorted (and indexed) BAM file')
    group4 = parser.add_argument_group("Bias calculation information","Highly recommended. If fasta is not provided, will not calculate bias")
    group4.add_argument('--fasta', metavar = 'genome_seq',
                    help = 'Indexed fasta file')
    group4.add_argument('--pwm', metavar = 'Tn5_PWM', help = "PWM descriptor file. Default is Human.PWM.txt included in package", default = "Human")
    group2 = parser.add_argument_group("General options","optional")
    group2.add_argument('--out', metavar = 'out_basename', help = "output file basename")
    group2.add_argument('--cores', metavar = 'num_cores',default=1,
                    help='Number of cores to use',type=int)
    group5 = parser.add_argument_group("NFR determination parameters")
    group5.add_argument('--max_occ', metavar= 'float', default = 0.1,
        help = 'Maximum mean occupancy for NFR. Default is 0.1', type = float)
    group5.add_argument('--max_occ_upper', metavar= 'float', default = 0.25,
        help = 'Maximum for minimum of  upper bound occupancy in NFR. Default is 0.25', type = float)
 

def add_vprocess_parser( subparsers):
    """Add argument parsers for the vprocess utility

    """
    parser = subparsers.add_parser("vprocess", help = "nucleoatac function:  Make processed vplot to use for nucleosome calling")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--out', metavar='output_basename',required=True)
    group2 = parser.add_argument_group('VPlot and Insert Size Options', 'Optional')
    group2.add_argument('--sizes', metavar='file' , help = 'Insert distribution file')
    group2.add_argument('--vplot', metavar='vmat_file',
                    help = 'Accepts VMat file.  Default is Vplot from S. Cer.',
                    default = nucleoatac.Magic.default_vplot)
    group3 = parser.add_argument_group('Size parameers', 'Use sensible values')
    group3.add_argument('--lower',metavar="int",default=105,
                    help="lower limit (inclusive) in insert size. default is 105",type=int)
    group3.add_argument('--upper',metavar="int",default=251,
                    help="upper limit (exclusive) in insert size. default 251",type=int)
    group3.add_argument('--flank',metavar="int",default=60,
                    help="distance on each side of dyad to include",type=int)
    group4 = parser.add_argument_group('Options', '')
    group4.add_argument('--smooth', metavar = "float", default = 0.75, type = float,
                    help="SD to use for gaussian smoothing.  Use 0 for no smoothing.")
    group4.add_argument('--plot_extra', action='store_true',default=False,
                    help="Make some additional plots")
    return

def add_nuc_parser( subparsers):
    """Add argument parsers for the nuc utility

    """
    parser = subparsers.add_parser("nuc", help = "nucleoatac function:  Call nucleosome positions and make signal tracks")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bed', metavar='bed_file' , help = 'Regions for which \
        to do stuff.', required=True)
    group1.add_argument('--vmat', metavar='vdensity_file', help = "VMat object", required=True)
    group1.add_argument('--bam', metavar='bam_file',
                    help = 'Accepts sorted BAM file', required=True)
    group1.add_argument('--out', metavar='basename',
                    help="give output basename", required = True)
    group2 = parser.add_argument_group('Bias options',"If --fasta not provided, bias not calculated")
    group2.add_argument('--fasta', metavar = 'genome_seq',
                    help = 'Indexed fasta file')
    group2.add_argument('--pwm', metavar = 'Tn5_PWM', help = "PWM descriptor file. Default is Human.PWM.txt included in package", default = "Human")
    group3 = parser.add_argument_group('General options', '')
    group3.add_argument('--sizes', metavar = 'fragmentsizes_file',
                       help = "File with fragment size distribution.  Use if don't want calculation of fragment size")
    group3.add_argument('--occ_track', metavar = 'occ_file',
                       help = "bgzip compressed, tabix-indexed bedgraph file with occcupancy track. Otherwise occ not determined for nuc positions.")
    group3.add_argument('--cores', metavar = 'num_cores',default=1,
                    help='Number of cores to use',type=int)
    group3.add_argument('--write_all', action="store_true", default = False,
                    help="write all tracks")
    group3.add_argument('--not_atac', dest = "atac",  action="store_false", default = True,
                    help="data is not atac-seq")
    group5 = parser.add_argument_group('Nucleosome calling parameters','Change with caution')
    group5.add_argument('--min_z', metavar='float', default = 3,
        help = 'Z-score threshold for nucleosome calls. Default is 3', type = float)
    group5.add_argument('--min_lr', metavar='float', default = 0,
        help = 'Log likelihood ratio threshold for nucleosome calls. Default is 0', type = float)
    group5.add_argument('--nuc_sep',metavar='int',default=120,type=int,
        help = "Minimum separation between non-redundant nucleosomes. Default is 120")
    group5.add_argument('--redundant_sep',metavar='int',default=25,type=int,
        help = "Minimum separation between redundant nucleosomes. Not recommended to be below 15. Default is 25")
    group5.add_argument('--sd',metavar='int', type=int,
        default=10, help = "Standard deviation for smoothing. Affect the \
        resolution at which nucleosomes can be positioned. Not recommended to \
        exceed 25 or to be smaller than 10. Default is 10" )
    return


def add_run_parser( subparsers):
    """Add argument parsers for the run utility

    """
    parser = subparsers.add_parser("run", help = "Main nucleoatac utility-- runs through occupancy determination & calling nuc positions")
    group1 = parser.add_argument_group('Required', 'Necessary arguments')
    group1.add_argument('--bed', metavar='bed_file' , help = 'Regions for which \
        to do stuff.', required=True)
    group1.add_argument('--bam', metavar='bam_file',
                    help = 'Accepts sorted BAM file', required=True)
    group1.add_argument('--out', metavar='output_basename',
                    help="give output basename", required=True)
    group1.add_argument('--fasta', metavar = 'genome_seq',
                    help = 'Indexed fasta file', required=True)
    group4 = parser.add_argument_group("Bias calculation parameters","")
    group4.add_argument('--pwm', metavar = 'Tn5_PWM', help = "PWM descriptor file. Default is Human.PWM.txt included in package", default = "Human")
    group3 = parser.add_argument_group('General options', '')
    group3.add_argument('--cores', metavar = 'num_cores',default=1,
                    help='Number of cores to use',type=int)
    group3.add_argument('--write_all', action="store_true", default = False,
                    help="write all tracks")
    return


