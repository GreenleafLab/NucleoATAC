# NucleoATAC
Python package for calling nucleosomes using ATAC-Seq data.

Non-Python Dependencies:  

* bedGraphToBigWig (from USCS Genome Browser:  http://hgdownload.cse.ucsc.edu/admin/exe/)
(must be on PATH)

Python Dependencies:

* python 2.7
* bx-python
* scipy
* numpy
* pysam
* cython
* matplotlib

Tested only on Ubuntu OS.  

Installation:

1. Clone repository.
2. Run `python setup.py install` within main directory
(use --user flag for user specific install.  See https://docs.python.org/2/install/ for more details)

Necessary files prior to calling nucleosomes:
* BAM file with aligned data.  Generally this will be filtered for high mapping quality, etc.
* Bed file with open chromatin regions.  This file must be sorted, and must not contain any overlapping regions!
* Tn5 bias track in bigWig format (or alternatively, gDNA insertions).  This must span all regions in bed file (with several hundred additional flanking bases)!


How to call nucleosomes with mostly default parameters:

```
nucleoatac run --bed <bed_file> --bam <bam_file> --bias <Tn5_bias_track> --out <output_basename> --cores <num_cores_to_use>
```

This command is actually running 3 different commands.  For greater versatility, run each separately.  To see options, use help:  

```
nucleoatac occ --help
```

```
nucleoatac vprocess --help
```
```
nucleoatac nuc --help
```


There is an additional nucleoatac call that can be made (to make a V-plot for specific sites):
```
nucleoatac vplot --help
```

