# NucleoATAC
Python package for calling nucleosomes using ATAC-Seq data.

Versions:  

* version 0 represents code as used for initial characterization of NucleoATAC method 
(e.g. as descriped in biorxiv manuscript).  This version may be a bit problematic in terms of the installation via pip or setup.py.  Adding main folder to $PYTHONPATH and bin folder to $PATH should allow scripts to be run.  Or 'pip install -e .' For this version checkout branch v0
* version 0.1+ involved extensive code reorganization, changes to input options, 
changes to way files are read or written, but no changes to algorithms.  Additional 
functions (not necessarily tied to nucleosome calling but useful for working with atac-seq data,
 are included under pyatac command line function).  Still undergoing some testing!

Instruction for use can be found at http://greenleaflab.github.io/NucleoATAC/

OS requirements:
Tested only on Ubuntu OS.

Python Dependencies:

* python 2.7
* scipy
* numpy
* pysam
* matplotlib
* cython

Installation:

Use pip -- `pip install .` inside NucleoATAC directory. 




