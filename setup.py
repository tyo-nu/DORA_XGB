#!/usr/bin/env python
from setuptools import setup,find_packages

def parse_requirements(filename):
    with open(filename, 'r') as file:
        return file.read().splitlines()

setup(
    name='DORA_XGB',

    version='1.0',

    description='Reaction feasibility models that take an input reaction string of the form A + B = C + D or of the form'
                'A.B>>C.D and output a feasibility score on the scale of 0 to 1 or a feasibility label (either 0 or 1).'
                'The DORA-XGB feasibility models are based on the gradient boosted trees (XGBoost) architecture.'
                'All hyperparameters were optimized based on a Bayesian hyperparameter optimization procedure.'
                'Our reaction feasibility models make predictions on the feasibility of a given reaction on two factors.'
                'The first is thermodynamic feasibility, which the model should have learned implicitly.'
                'And the second is the likelihood of an alternate reaction center undergoing an enzymatic reaction.',

    author = "Yash Chainani and Joseph Ni",

    author_email = "yashchainani2026@u.northwestern.edu",

    packages=find_packages(),

    install_requires=parse_requirements('requirements.txt') + ['map4 @ git+https://github.com/reymond-group/map4@v1.0'],

    package_dir = {},

    package_data = {'feasibility_wrapper': ['models/*',
                                            'cofactors/expanded_cofactors_no_stereochem.tsv']},

classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research/Development",
        "Intended Audience :: Scientific Engineering",
        "Intended Audience :: Application",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
    ],
)