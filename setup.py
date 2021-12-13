#!/usr/bin/env python3

import setuptools
import os

from distutils.extension import Extension

if __name__ == "__main__":
    setuptools.setup(
        name="panorama",
        version=open(os.path.join(os.path.dirname(__file__), "VERSION")).read().rstrip(),
        # url="https://github.com/labgem/",
        description="Comparative pangenomic analysis toolbox",
        packages=setuptools.find_packages(),
        setup_requires=["cython"],
        install_requires=[],
        package_data={'': ['rRNA_DB/*cm*']},
        classifiers=["Environment :: Console",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
                     "Natural Language :: English",
                     "Operating System :: POSIX :: Linux",
                     "Programming Language :: Python :: 3",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"],
        entry_points={"console_scripts": ["panorama = panorama.main:main"]},
    )
