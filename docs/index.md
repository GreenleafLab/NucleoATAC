# NucleoATAC Documentation

For code visit [GitHub](https://github.com/GreenleafLab/NucleoATAC).

## About NucleoATAC

NucleoATAC is a python package for calling nucleosome positions and occupancy using ATAC-Seq data.  Functions for calling nucleosomes are included in the [nucleoatac](nucleoatac.md) command-line function.  NucleoATAC also includes other utlities for working with ATAC-seq data under the [pyatac](pyatac.md) function. 

## QuickStart

###Needed files
* Aligned paired-end reads in BAM format.  Must be sorted & indexed.  Probably should be filtered for quality.
* Fasta file with genome reference.  Must be indexed by faidx from [samtools](http://www.htslib.org/).
* Bed file with regions to perform analysis.  Generally will be broad open chromatin regions.

###Installation
Recommended to install NucleoATAC within virtual environment to minimize chance of issues with older versions of python package dependencies. NucleoATAC requries Python 2.7 (Does not currently support Python 3).

* Make and activate python [virtual environment](https://virtualenv.pypa.io/en/latest/)
* Use pip to install python dependencies: `pip install cython numpy scipy matplotlib pysam`

Either within a python virtual environment or system-wide:

* Get source code:  `git clone https://github.com/GreenleafLab/NucleoATAC`
* Change into directory: `cd /path/to/NucleoATAC`
* Install:  `pip install .`

For more installation options and troubleshooting, see [Installation page](installation.md)

###Calling nucleosomes
`nucleoatac run --bed <bedfile> --bam <bamfile> --fasta <fastafile> --out <output_basename> --cores <number_cores_to_use>`

Run `nucleoatac run --help` for all options.  For outputs and more details see [nucleoatac page](nucleoatac.md)



