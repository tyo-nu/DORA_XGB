#!/usr/bin/env python
from setuptools import setup,find_packages

def parse_requirements(filename):
    with open(filename, 'r') as file:
        return file.read().splitlines()

setup(
    name='DORA_XGB',

    version='1.10',

    description = "Gradient-booseted classifiers to predict the feasibility of enzymatic reactions",
                
    author = "Yash Chainani and Joseph Ni",

    author_email = "yashchainani2026@u.northwestern.edu",

    packages=find_packages(),

    install_requires=parse_requirements('requirements.txt'),

    package_dir = {},

    include_package_data = True,

    package_data = {'DORA_XGB': ['models/*',
                                 'cofactors/expanded_cofactors_no_stereochem.tsv']},

classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
    ],
)
