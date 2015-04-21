# NucleoATAC Documentation

For code visit [GitHub](https://github.com/GreenleafLab/NucleoATAC).

## About NucleoATAC

NucleoATAC is a python package for calling nucleosome positions and occupancy using ATAC-Seq data.  NucleoATAC also includes other utlities for working with ATAC-seq data.  

## QuickStart

###Needed files
* Aligned paired-end reads in BAM format.  Must be sorted & indexed.  Probably should be filtered for quality.
* Fasta file with genome reference.  Must be indexed by faidx from samtools.
* Bed file with regions to perform analysis.  Generally will be broad open chromatin regions.

###Installation
* Get source code:  `git clone https://github.com/GreenleafLab/NucleoATAC`
* Change into directory: `cd /path/to/NucleoATAC`
* Install:  `pip install .`

###Calling nucleosomes
`nucleoatac run --bed <bedfile> --bam <bamfile> --fasta <fastafile> --out <output_basename> --cores <number_cores_to_use>`

Run `nucleoatac run --help` for all options.  



