#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import setup
import os
import pysplicer

with open(os.path.join(os.path.dirname(__file__), "README.md")) as readme:
    long_description = readme.read()

packages = [
    "pysplicer"
]

setup(
    name="pysplicer",
    version=pysplicer.__version__,
    url="https://gitorious.org/pysplicer/pysplicer",
    license="GNU Affero General Public License v3",
    author="Cathal Garvey",
    author_email="cathalgarvey@cathalgarvey.me",
    maintainer="Cathal Garvey",
    maintainer_email="cathalgarvey@cathalgarvey.me",
    description="A frequency-bias codon optimisation script with early NGG avoidance and optional IUPAC match avoidance.",
    long_description=long_description,
    packages=packages,
    install_requires=[],
    scripts=["bin/pysplicer", "bin/cud-to-pysplicer"],
    platforms="any",
    zip_safe=False,
    classifiers=[
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Development Status :: 3 - Alpha",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
        "Environment :: Console",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.1",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
    ]
)
