#!/usr/bin/env python

##### IMPORT MODULES #####
# import necessary for python
from __future__ import division
import os
import sys
import string
import datetime
from optparse import OptionParser
import numpy as np
import gzip
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement

#### OPTIONS ####
# read options from command line
opts = OptionParser()
opts.add_option("-i", "--input", default="", 
help="<inputFile> Accepts FASTA file with single Sequence")
opts.add_option("-p", "--pwm", default="",
help="<pwmFile> File with normalized Tn5 transposon preference. Insertion point should be 11th position")
opts.add_option("-c", "--chrN", default="",
help="<chrN> name of chromosome (i.e. chr1)")
opts.add_option("-o", "--out", default="",
help="<outFile> Base name for Output file")
options, arguments = opts.parse_args()

##### INPUTS AND OUTPUTS #####
in_file = options.input
pwm_file = options.pwm
chrN=options.chrN
a=10#bases around transposition center included in motif

if options.input == "":
    sys.exit("Requires an input file, use --input to designate. Accepts FASTA files.")

if options.out=="":
	out = os.path.splitext(os.path.basename(in_file))[0]
else:
	out=options.out


#open PWM file
pwm=np.loadtxt(pwm_file)
w=len(pwm[0])
pwm=np.log(pwm)

##write function to determine the pwm score...
def scorer_mono(seq,m):
	#m should be array of size 4 * len(seq)
	S=0
	for i in range(len(seq)):
		if seq[i]=='A':
			S=S+m[0][i]
		elif seq[i]=='C':
			S=S+m[1][i]
		elif seq[i]=='G':
			S=S+m[2][i]
		elif seq[i]=='T':
			S=S+m[3][i]
	#logS=np.log(S)
	return S

def scorer_di(seq,m):
	#m should be array of size 16 * len(seq)
	S=0
	for i in range((len(seq)-1)):
		if seq[i:(i+2)]=='AA':
			S=S+m[0][i]
		elif seq[i:(i+2)]=='AC':
			S=S+m[1][i]
		elif seq[i:(i+2)]=='AG':
			S=S+m[2][i]
		elif seq[i:(i+2)]=='AT':
			S=S+m[3][i]
		elif seq[i:(i+2)]=='CA':
			S=S+m[4][i]
		elif seq[i:(i+2)]=='CC':
			S=S+m[5][i]
		elif seq[i:(i+2)]=='CG':
			S=S+m[6][i]
		elif seq[i:(i+2)]=='CT':
			S=S+m[7][i]
		elif seq[i:(i+2)]=='GA':
			S=S+m[8][i]
		elif seq[i:(i+2)]=='GC':
			S=S+m[9][i]
		elif seq[i:(i+2)]=='GG':
			S=S+m[10][i]
		elif seq[i:(i+2)]=='GT':
			S=S+m[11][i]
		elif seq[i:(i+2)]=='TA':
			S=S+m[12][i]
		elif seq[i:(i+2)]=='TC':
			S=S+m[13][i]
		elif seq[i:(i+2)]=='TG':
			S=S+m[14][i]
		elif seq[i:(i+2)]=='TT':
			S=S+m[15][i]
	#logS=np.log(S)
	return S

def Scorer(seq,m):
	if len(m)==4:
		return scorer_mono(seq,m)
	elif len(m)==16:
		return scorer_di(seq,m)

#Get sequence to score 

#if no chromosome is given, just use first sequence in fasta input
if chrN=='':
	rec=SeqIO.parse(in_file,'fasta').next()
	source=rec.seq
	chrN=rec.id
else:
	for rec in SeqIO.parse(in_file,'fasta'):
		if rec.id==chrN:
			source=rec.seq

Seq_length=len(source)

#open output file for writing scores
if options.chrN!='':
	output=open(out+'.'+chrN+'.Scores.bed','w')
else:
	output=open(out+'.Scores.bed','w')


for i in range(a,Seq_length-w+a):
	forwardseq=str(source[i-a:i+w-a]).upper()
	f=Scorer(forwardseq,pwm)
	output.write('\t'.join([chrN,str(i),str(i+1),str(f)])+'\n')


output.close()


