# MDSubSampler: Molecular Dynamics SubSampler

[![PyPI version](https://badge.fury.io/py/mdsubsampler.svg)](https://badge.fury.io/py/mdsubsampler)

MDSubSampler is a Python library and toolkit for a posteriori subsampling of multiple trajectory data for further analysis. This toolkit implements uniform, random, stratified sampling, bootstrapping and targeted sampling to preserve the original distribution of relevant geometrical properties.

## Prerequisites

This project requires Python (version 3.9.1 or later). To make sure you have the right version available on your machine, try running the following command. 

```sh
$ python --version
Python 3.9.1
```

## Table of contents

- [Project Name](#project-name)
  - [Prerequisites](#prerequisites)
  - [Table of contents](#table-of-contents)
  - [Getting Started](#getting-started)
  - [Installation](#installation)
  - [Usage](#usage)
    - [Workflow](#workflow)
    - [Scenarios](#scenarios)
    - [Parser](#parser)
    - [Development](#development)
  - [Authors](#authors)
  - [License](#license)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for analysis and development purposes. 

## Installation

**BEFORE YOU INSTALL:** please read the [prerequisites](#prerequisites)

To install and set up the library, run:

```sh
pip install MDSubSampler
```

## Usage 

### Workflow

Input:
- Molecular Dynamics trajectory 
- Geometric property
- Atom selection [optional - default is "name CA"]
- Reference structure [optional] 
- Sample size or range of sizes
- Dissimilarity measure [optional - default is "Bhattacharyya"]

Output:
- .dat file with calculated property for full trajectory (user input)
- .dat file(s) with calculated property for one or all sample sizes input
- .xtc file(s) with sample trajectory for one or all sample sizes
- .npy file(s) with sample trajectory for one or all sample sizes 
- .npy training set for ML purposes for sample trajectory (optional)
- .npy testing set for ML purposes for sample trajectory (optional)
- .npy file(s) with sample trajectory for one or for all sample sizes 
- .png file with overlapped property distribution of reference and sample
- .json file report with important statistics from the analysis
- .txt log file with essential analysis steps and information

### Scenarios

To run scenarios 1,2 or 3 you can download your protein trajectory and topology file (.xtc and .gro files) to the data folder and then run the following:

```sh
python mdss/scenarios/scenario_1.py data/<YourTrajectoryFile>.xtc data/<YourTopologyfile>.gro <YourPrefix>
```
Scenarios 1,2 and 3 are also available in Jupyter Notebooks format, can be used as templates and can be modified interactively according to the user's needs. You can also find more advanced scenarios in the cookbook directory. If you clone the library locally to your machine, then run the following command before you run the cells.  

```sh
%cd <pathToMDSubSamplerDirectory>
```

### Parser

If you are a terminal lover you can use the terminal to run the code and make a choice for the parser arguments. To see all options and choices run:

```sh
python mdss/run.py --help
```
Once you have made a selection of arguments, your command can look like the following example:

```sh
python mdss/run.py \
  --traj "data/<YourTrajectoryFile>.xtc" \
  --top "data/<YourTopologyFile>.gro" \
  --prefix "<YourPrefix>" \
  --output-folder "data/<YourResultsFolder>" \
  --property='DistanceBetweenAtoms' \
  --atom-selection='G55,P127' \
  --sampler='BootstrappingSampler' \
  --n-iterations=50 \
  --size=<SampleSize> \
  --dissimilarity='Bhattacharyya'
```

### Development

#### With Poetry

Start by either downloading the tarball file from https://github.com/alepandini/MDSubSampler to your local machine or cloning this repo on your local machine:

```sh
git clone git@github.com:alepandini/MDSubSampler.git
cd MDSubSampler
```

Following that, download and install poetry from https://python-poetry.org/docs/#installation


Finally, run the following:

```sh
poetry install
poetry build
poetry shell
```
You can now start developing the library.

#### With Docker

Start by installing Docker using this link https://docs.docker.com/get-docker/.

Initially a Docker image will need to be built. To do this run the following command:

```sh
docker build -t <image name> .
```

Then run the following command to get access to a shell with all dependencies installed:

```sh
docker run -it -v $(pwd):/app -e PYTHONPATH=/app <image name> /bin/bash
```

This will also mirror the local filesystem in the Docker image, so that any local change will
be reflected in the running container, and vice-versa, using a Docker volume.

The repo also includes two handy scripts to run all of the above faster (an image called
`subsampler` will be created):

```sh
./build-docker
./run-docker
```

After dropping in the Docker shell, all dependencies will be installed, and the package scripts
will also be in scope (the `mdss` command and all scenarios declared in `pyproject.toml`).

### Authors

* **Namir Oues** - [namiroues](https://github.com/namiroues)
* **Alessandro Pandini** [alepandini](https://github.com/alepandini)

### License

The library is licensed by **GPL-3.0**
