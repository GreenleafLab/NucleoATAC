from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys
import os

class NoTestCommand(TestCommand):
    def run(self):
        print("NucleoATAC does not support running tests with "
              "'python setup.py test'. The test suite can be run via 'pytest'")


setup(name='NucleoATAC',
    version='0.4.0',
    description='python package for calling nucleosomes with ATAC-Seq',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='ATAC-Seq sequencing bioinformatics',
    url='https://github.com/GreenleafLab/NucleoATAC',
    author='Alicia Schep',
    author_email='aschep@gmail.com',
    license='MIT',
    packages=['pyatac','pyatac.pwm','nucleoatac','nucleoatac.vplot'],
    python_requires='>=3.8',
    install_requires=['cython > 0.22', 'numpy >= 1.9.1', 'scipy >= 0.16.0','pysam >= 0.10.0','matplotlib'],
    scripts=['bin/pyatac','bin/nucleoatac'],
    include_package_data=True,
    zip_safe=False,
    cmdclass = {'test': NoTestCommand})
