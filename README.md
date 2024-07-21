# DORA-XGB: An improved enzymatic reaction feasibility classifier trained using a novel synthetic data approach

Authors: Yash Chainani, Zhuofu Ni, Kevin M. Shebek, Linda J. Broadbelt, and Keith E.J. Tyo

This public repository holds the supervised learning reaction feasibility models, DORA_XGB. For examples on how to use our models, see `scripts/run_example.py` or  `notebooks/DORA_XGB_examples.ipynb`.

## Environment setup and installation:
To use DORA-XGB models, begin by creating a python 3.8 environment or use an existing python 3.8 environment:  

```
conda create -n DORA_XGB_env python=3.8
```

After creating a new python 3.8 virtual environment, or using an existing one, simply pip install our DORA-XGB models package:  

```
pip install -e.
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
