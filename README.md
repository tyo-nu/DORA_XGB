# DORA-XGB: An improved enzymatic reaction feasibility classifier trained using a novel synthetic data approach

[Chainani, Yash, Zhuofu Ni, Kevin M. Shebek, Linda J. Broadbelt, and Keith EJ Tyo. "DORA-XGB: an improved enzymatic reaction feasibility classifier trained using a novel synthetic data approach." Molecular Systems Design & Engineering (2025)](https://pubs.rsc.org/en/content/articlehtml/2024/me/d4me00118d)

This public repository holds the supervised learning reaction feasibility models, DORA_XGB. For examples on how to use our models, see `scripts/run_example.py` or  `notebooks/DORA_XGB_examples.ipynb`.

## Environment setup & installation options:
#### 1. Clone this repository:
To use our DORA-XGB models, clone this repository then create a new python environment or use an existing one. Subsequently, from the same directory as the `setup.py` file, pip install DORA-XGB:
```
conda create -n DORA_XGB_env python=3.8
pip install -e.
```
#### 2. Install from the PyPI repository:
The most convenient way to begin using our DORA-XGB models may be to directly install them from the python package index (PyPI):
```
pip install DORA-XGB
```

## Running DORA-XGB with docker:
We have also created a docker container for users to deploy our models within a containerized environment. To begin, run the following in the same directory as the dockerfile to build a docker image with the name `dora_xgb`:  

`docker build -t dora_xgb .`

After building the docker image locally, spin up a container with an interactive bash shell:

`docker run -ti dora_xgb /bin/bash`

In this interactive bash shell, the run_example script and be run simply using:

`python run_example.py`

To edit the contents of each script, you can download the vim text editor in the docker container:

`apt-get update && apt-get install -y vim`

To shut down a docker container and return to your terminal, simply type `exit` into the interactive bash shell.
