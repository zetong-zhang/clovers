"""
Setup script for building the Zcurve Python extension module.

This module provides DNA sequence encoding functionality using the Z-curve method
for bioinformatics applications, particularly useful in gene prediction.
"""
from setuptools import setup, Extension
import numpy as np


# Extension module definition
zcurve_extension = Extension(
    "Zcurve",
    sources=["Zcurve.cpp"],
    include_dirs=[
        ".",
        np.get_include(),
    ],
)

setup(
    name="Zcurve",
    version="0.0.1",
    description="DNA sequence encoding using Z-curve method",
    long_description=__doc__,
    author="Zetong Zhang",
    author_email="zhangzetong@tju.edu.cn",
    url="https://github.com/zetong-zhang/clovers",
    license="GPLv3",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
    ],
    keywords="bioinformatics genomics z-curve gene-prediction dna-sequence",
    packages=[],
    ext_modules=[zcurve_extension],
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=["numpy>=1.26.0"]
)
