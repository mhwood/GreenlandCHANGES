# Greenland CHANGES
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

Next, install the dependencies of the GreenlandCHANGES package as
```
conda install -c conda-forge gdal
pip install shapely
conda install -c anaconda geos
conda install -c conda-forge scipy
pip install numpy
pip install netcdf4
pip install xarray
pip install pyproj==2.0.0
```

Finally, the GreenlandCHANGES package can be cloned as 
```
git clone https://github.com/mhwood/GreenlandCHANGES.git
```

## Using the GreenlandCHANGES package
For convenience, we provide several tutorials, presented as jupyter notebooks, in the `tutorials` directory.
