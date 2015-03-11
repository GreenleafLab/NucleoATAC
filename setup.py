from setuptools import setup
import sys


def readme():
    with open('README.md') as f:
        return f.read()

if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)

setup(name='nucleoatac',
    version='0',
    description='nucleosome calling for ATAC-Seq',
    long_description=readme(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='nucleosome ATAC-Seq sequencing',
    url='https://github.com/GreenleafLab/NucleoATAC',
    author='Alicia Schep',
    author_email='aschep@stanford.edu',
    license='MIT',
    packages=['nucleoatac'],
    install_requires=['numpy', 'scipy','pysam','matplotlib','cython','bx-python'],
    scripts=['bin/nucleoatac'],
    include_package_data=True,
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'])



