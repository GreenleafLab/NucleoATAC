###Dependencies 
NucleoATAC was developed and tested on Ubuntu (14.04.1 LTS).
Not designed for compatability with Windows.

Python modules required:

* numpy (>=1.9.1)
* scipy
* matplotlib
* cython (>=0.22)
* pysam (>= 0.8.1)


###Installation

1) Either download latest release or clone repository to ensure you get the latest version: 

*  See [releases](https://github.com/GreenleafLab/NucleoATAC/releases)
* `git clone https://github.com/GreenleafLab/NucleoATAC.git`

2) Either:

* Use [pip](https://pip.pypa.io/en/latest/) to install, i.e. `pip install .` when in the NucleoATAC directory
* Use [setup.py](https://docs.python.org/2/install/), i.e. `python install setup.py` when in the NucleoATAC directory

For both, installation can be user-specific if desired.  See documentation for pip or setup.py for options on how to do user-specific installation or other features of either installation method.

###Virtual Environment

Installing nucleoatac within a virtual environment is recommended!  Use [virtualenv](https://virtualenv.pypa.io/en/latest/) to create a virtual environment.  For some reason I've found that installing some of the dependencies (mainly scipy) prior to installing NucleoATAC (which will try to install all dependencies) seems to work best.  

###Troubleshooting

What is version of cython?  If you run `import cython` then `cython.__version__` is version 0.22?  I've found that even with the requirement for cython >= 0.22 in setup.py that NucleoATAC may install even if python doesn't automatically import cython v0.22.

Try installing in virtual environment (see above). Sometimes if older versions of some python packages are available it can be difficult to cleanly upgrade; using a virtual environment allows you to freshly install necessary packages.



