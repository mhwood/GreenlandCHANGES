# GreenlandCHANGES
**C**omprehensive **H**omogenization **A**nd dow**N**loading of **G**lacier **E**levation and **S**peed data

![Jakobshavn Example](https://github.com/mhwood/GreenlandCHANGES/blob/master/jakobshavn_example.jpg)

## Motivation
Today, there is an ever-growing volume of ice elevation and velocity measurements on glaciers around the Greenland Ice Sheet. These unprecedented data sources represent a key cornerstone of scientific research into how the ice sheet is losing mass and contributing to sea level rise. However, identifying pertinent data and preparing it for analysis represents a time-consuming, yet routine task standing in the way of an initial research idea and tangible results. **The goal of this package is to provide researchers with a convenient way to download and homogenize all available ice velocity and/or elevation data in regions of interest in a format which is convenient for analysis.**

## Overview of the GreelandCHANGES package

The main functionality of this package rests in the core `changes` module. This module identifies, downloads, homogenizes, and compiles all available ice velocity and elevation data in a user-defined region of interest. The output product is a list of data "stacks" - a compilation of homogenized data stored on identical grids with the same geographic projection and resolution. A concise description of this module is available in the changes module [README](https://github.com/mhwood/GreenlandCHANGES/blob/master/changes/README.md).

This package also provides several tools that facilitate simple analysis and plotting of data produced by the `changes` module. These tools are available in the `toolbox`.

## :exclamation: Important User Warning
The data sets utilized in this package are from a variety of sources, each of which must cited appropriately when they are used in publications, conference presentations, or in any other form. We have provided the listed citation for each source in the README files for both [elevation](https://github.com/mhwood/GreenlandCHANGES/blob/master/changes/elevation/README.md) and ice [velocity](https://github.com/mhwood/GreenlandCHANGES/blob/master/changes/velocity/README.md). 

## Getting started
It is recommended that a fresh anaconda environment is used for this code. For example, a new environment called `greenlandchanges` can be created and activated as
```
conda create --name greenlandchanges
conda activate greenlandchanges
```

Next, install the dependencies of the GreenlandCHANGES package and the tutorials herein:
```
conda install -c anaconda netcdf4
conda install -c conda-forge matplotlib
conda install -c anaconda xarray
conda install -c conda-forge pyproj
conda install -c conda-forge scipy
conda install -c conda-forge jupyterlab
conda install -c anaconda ipykernel
python3 -m ipykernel install --user --name=greenlandchanges
conda install -c conda-forge gdal
```
This order checked on MacOS 12 Oct 21.

Finally, the GreenlandCHANGES package can be cloned as 
```
git clone https://github.com/mhwood/GreenlandCHANGES.git
```

## Using the GreenlandCHANGES package
For convenience, we provide several tutorials, presented as jupyter notebooks, in the [tutorials](https://github.com/mhwood/GreenlandCHANGES/tree/master/tutorials) directory.

To use your conda environment in the jupyter notebook tutorials, make sure jupyter is installed in your environment
```
conda install -c conda-forge jupyterlab
```
and ensure that it is available in the IPy kernel:
```
conda install -c anaconda ipykernel
python3 -m ipykernel install --user --name=[name of environment]
```

## How To Configure Earthdata Login Authentication
Follow steps 1-3 to configure access to EarthData:

[Step-by-step Guide to Configure Authentication Credentials](https://urs.earthdata.nasa.gov/documentation/for_users/data_access/curl_and_wget)

