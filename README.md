# MDSubSampler: Molecular Dynamics SubSampler

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
  - [Contributing](#contributing)
  - [Credits](#credits)
  - [Built With](#built-with)
  - [Versioning](#versioning)
  - [Authors](#authors)
  - [License](#license)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for analysis and development purposes. 

## Installation

**BEFORE YOU INSTALL:** please read the [prerequisites](#prerequisites)

Start by either downloading the turball file from https://github.com/alepandini/MDSubSampler to your local machine or cloning this repo on your local machine:

```sh
$ git clone git@github.com:alepandini/MDSubSampler.git
$ cd MDSubSampler
```

To install and set up the library, run:

```sh
$ pip install MDSubSampler
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

To run scenarios 1,2 or 3 you can download your protein trajectory to the data folder and then run the following:

```sh
$ python mdss/scenarios/scenario_1.py data/<YourInputXTCFile> data/<YourInputPDBFile> <YourPrefix>
```

### Parser

If you are a terminal lover you can use the terminal to run the code and make a choice for the parser arguments. For example you can run the following:


```sh
$ python mdss/run.py --traj "data/input.xtc" --top "data/input.gro" --prefix "<YourPrefix>" --output-folder "data/results" --property='DistanceBetweenAtoms' --atom-selection='G55,P127' --sampler='BootstrappingSampler' --n-iterations=50 --size=11000 --dissimilarity='Bhattacharyya'
```

### Development

Finally, if you want to develop the library further you can start by installing poetry:

```sh
$ curl -sSL https://install.python-poetry.org | python3 -
```

Then you can run the following:

```sh
$ poetry install
Installing dependencies from lock file
No dependencies to install or update
Installing the current project: mdsubsampler (0.1.0)
$ poetry build
Building mdsubsampler (0.1.0)
  - Building sdist
  - Built mdsubsampler-0.1.0.tar.gz
  - Building wheel
  - Built mdsubsampler-0.1.0-py3-none-any.whl
$ poetry shell
(mdsubsampler-py3.9) ~/
```
Following that you can start developing the code and commit changes to GitHub. 

### Authors

* **Namir Oues** - [namiroues](https://github.com/namiroues)
* **Alessandro Pandini** [alepandini](https://github.com/alepandini)


