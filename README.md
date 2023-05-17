# aetherpy
[![Coverage Status](https://coveralls.io/repos/github/AetherModel/aetherpy/badge.svg?branch=main)](https://coveralls.io/github/AetherModel/aetherpy?branch=main)

The Python package that supports Aether model data management and analysis.

# Installation

## Starting from scratch
* Python and the packages aetherpy depends on are freely available. Python may
  be obtained through various package managers or from https://www.python.org/.
  A common OS-independent Python package manager to obtain the aetherpy
  dependencies is [PyPi](https://pypi.org/).

## Installation from GitHub

```
git clone https://github.com/aaronjridley/aetherpy
cd aetherpy
git checkout develop

python setup.py install --user
# <or>
python setup.py develop --user
```

If you want to install aetherpy for the entire system, you will need
to take off the `--user` flag and add `sudo` in front of the
command in the above block:
```
sudo python setup.py install
```

The difference between the `install` and `develop` methods is that
develop allows you to change branches to the desired test brach or
alter files without needing to re-install aetherpy.

## Installation from PyPi

Pip installation will soon be available

# Getting Started

aetherpy contains a test script to plot model results.  It may be called via:

```
python aetherpy/run_plot_model_results.py var=3 -alt=120 [file_dir]/3DALL*.bin
```

where [file_dir] is the directory where you have Aether model output files. You
can see the help by typing:

```
python aetherpy/run_plot_model_results.py -h
```
