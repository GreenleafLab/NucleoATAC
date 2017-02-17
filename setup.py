from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys
import os

if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)


class NoTestCommand(TestCommand):
    def run(self):
        print("NucleoATAC does not support running tests with "
              "'python setup.py test'. Please run 'python tests.py'")


setup(name='NucleoATAC',
    version='0.3.4',
    description='python package for calling nucleosomes with ATAC-Seq',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='ATAC-Seq sequencing bioinformatics',
    url='https://github.com/GreenleafLab/NucleoATAC',
    author='Alicia Schep',
    author_email='aschep@stanford.edu',
    license='MIT',
    packages=['pyatac','pyatac.pwm','nucleoatac','nucleoatac.vplot'],
    install_requires=['cython >= 0.22','numpy >= 1.9.1', 'scipy >= 0.16.0','pysam >= 0.10.0','matplotlib'],
    scripts=['bin/pyatac','bin/nucleoatac'],
    include_package_data=True,
    zip_safe=False,
    tests_require=['nose'],
    cmdclass = {'test': NoTestCommand})
