from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "pyatac.fragments",
        ["pyatac/fragments.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "nucleoatac.multinomial_cov",
        ["nucleoatac/multinomial_cov.pyx"],
        include_dirs=[np.get_include()],
    ),
]

setup(
    name="NucleoATAC",
    version="1.0.0",
    description="Nucleosome calling from ATAC-seq",
    author="Alicia Schep",
    author_email="aschep@stanford.edu",
    url="https://github.com/GreenleafLab/NucleoATAC",
    license="MIT",

    packages=find_packages(),

    python_requires=">=3.11",

    install_requires=[
        "numpy>=1.26",
        "scipy>=1.11",
        "pysam>=0.22",
        "matplotlib>=3.8"
    ],

    ext_modules=cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
    ),

    include_dirs=[np.get_include()],

    entry_points={
        "console_scripts": [
            "nucleoatac=nucleoatac.cli:main",
            "pyatac=pyatac.cli:main",
        ],
    },

    include_package_data=True,
    zip_safe=False,

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
