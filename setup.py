'''
HPVMod installation. Requirements are listed in requirements.txt.
'''

import os
import sys
import runpy
from setuptools import setup, find_packages

# Load requirements from txt file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Get version
cwd = os.path.abspath(os.path.dirname(__file__))
versionpath = os.path.join(cwd, 'hpvmod', 'version.py')
version = runpy.run_path(versionpath)['__version__']

# Get the documentation
with open(os.path.join(cwd, 'README.rst'), "r") as f:
    long_description = f.read()

CLASSIFIERS = [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: Other/Proprietary License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Development Status :: 5 - Production/Stable",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
]

setup(
    name="hpvmod",
    version=version,
    author="Robyn Stuart and others",
    author_email="robyna.s@gmail.com",
    description="HPV Model",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url='',
    keywords=["HPV", "compartmental model", "interventions", "epidemiology"],
    platforms=["OS Independent"],
    classifiers=CLASSIFIERS,
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    extras_require={
        'full':  [
            'plotly',
            'fire',
            'optuna',
            'parestlib',
        ],
    }
)

