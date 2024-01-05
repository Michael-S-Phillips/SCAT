# SpectralCubeAnalysisTool

## Overview

The `SpectralCubeAnalysisTool` class is designed to provide a graphical user interface (GUI) for analyzing multi- or hyperspectral data. It allows users to load hyperspectral data and analyze it, including displaying RGB images, drawing and analyzing polygonal ROIs on the data (and exporting to GIS compatible shape layers), and extracting spectral plots from points and ROIs. This tool was originally designed for analysis of CRISM data, therefore spectral ratioing capabilities are built in. 

The intended workflow is to load the hyperspectral data into the left image display and to load a reduced dimensionality representation into the right image display (such as spectral parameters or PCA/MNF/ICA bands). This allows the user to select ROIs in the right image display and extract spectra from those same regions from the spectral cube.

This tool is a work in progress! I will be updating this tool as I work on it. Currenlty, I've only used this for testing purposes and have only loaded in a few different types of CRISM data. The tool only has functionality to read in ENVI .img files with a .hdr associated with them. Perhaps in the future I will expand this to read in a more diverse set of image classes.

## Features

### RGB Display

- The class allows users to visualize hyperspectral data in RGB format using user-defined bands.
- Users can adjust the color map stretch for the RGB display.

### Drawing and Analyzing Polygons

- The tool enables users to draw and manipulate polygon shape layers on the data.
- Users can extract spectral data statistics (average and standard deviation) from the drawn polygons.
- Geospatial information associated with the drawn polygons is also accessible.

### Interface Elements

- The GUI provides buttons for resetting the display extent, toggling the polygon drawing mode, clearing polygons from the display, redrawing existing polygons, deleting all polygons, and plotting mean spectra of polygons.
- Users can load hyperspectral cubes and band parameter images.
- Histogram plotting options are available for both the left and right frames.

## Getting Started

# Installaing
```bash
git clone git@github.com:Michael-S-Phillips/SCAT.git
conda env create -f environment.yml
```
## Dependencies

- The code relies on the `matplotlib`, `tkinter`, and `numpy` libraries for graphical and numerical operations. Other libraries are listed in the environment.yml file.

## How to Run

from the command line:
```bash
python scat.py
```
To use the `SpectralCubeAnalysisTool` class, you can follow these steps:

1. Load hyperspectral data and a band parameter image using the provided menu options (file>Load Hyperspectral Cube).
2. Adjust the display settings, such as band selection and stretch.
3. Use the buttons in the interface to draw polygons, save polygons as geo shape files, and analyze data.

## Incorporating HyPyRameter for Spectral Parameter Calculation

SCAT provides a GUI interface for using the HyPyRameter library for spectral parameter calculation. Go to the HyPyRameter page to learn more about this tool (https://github.com/Michael-S-Phillips/HyPyRameter). 


## TODO List

The code includes a TODO list of features and improvements to be implemented in the future:

- Adding more functionality to the polygon shape layer drawing.
- Spectral smoothing routines
- Adding basic spectral analysis workflows, such as MNF