from setuptools import setup
import sys


if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)


setup(name='NucleoATAC',
    version='0.2.1',
    description='python package for calling nucleosomes with ATAC-Seq',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='ATAC-Seq sequencing bioinformatics',
    url='https://github.com/GreenleafLab/PyAtac',
    author='Alicia Schep',
    author_email='aschep@stanford.edu',
    license='MIT',
    packages=['pyatac','pyatac.pwm','nucleoatac','nucleoatac.vplot'],
    install_requires=['cython >= 0.22','numpy', 'scipy','pysam >= 0.8.1','matplotlib'],
    scripts=['bin/pyatac','bin/nucleoatac'],
    include_package_data=True,
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'])



