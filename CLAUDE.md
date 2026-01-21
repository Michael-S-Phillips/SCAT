# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SCAT (Spectral Cube Analysis Tool) is a Python/Tkinter GUI application for analyzing hyperspectral imaging data. Originally designed for CRISM (Mars orbital spectrometer) data, it supports any ENVI-format hyperspectral imagery.

**Core workflow**: Load a hyperspectral cube in the left display, load a reduced-dimensionality representation (spectral parameters, PCA/MNF bands) in the right display, draw ROIs on either display, and extract/analyze spectra from the cube.

## Commands

```bash
# Setup
conda env create -f environment.yml
conda activate scat

# Run the application
python scat.py
```

## Architecture

The application is a monolithic Tkinter GUI in `scat.py` (~5500 lines) with the main class `SpectralCubeAnalysisTool`.

**Key components:**
- **Dual canvas display**: Left canvas shows full hyperspectral cube (RGB composite), right canvas shows reduced dimensionality data. Both support independent RGB band selection and histogram stretching.
- **ROI system**: Polygon drawing, statistics extraction (mean/std spectra), export to GeoPackage/Shapefile.
- **Spectral extraction**: Point-and-click spectral plots, optional ratioing (for CRISM mineral detection).
- **HyPyRameter integration**: GUI wrapper for spectral parameter calculation.

**Helper modules:**
- `create_histograms.py` - Interactive histogram windows for stretch adjustment
- `create_spectral_plot.py` - Spectral plot windows extending the main class

**External dependency**: `hypyrameter` library (from conda channel `michael--s--phillips`) for spectral parameter calculations.

## Data Formats

- **Input**: ENVI `.img` files with associated `.hdr` header files
- **Spectral libraries**: SPLIB format files in `librarySpectra/` directory (30+ mineral spectra)
- **ROI export**: GeoPackage, Shapefile (via GeoPandas)

## Development Guidelines

See `.claude/agents/` for specialized guidance on:
- `agent-spectral.md` - Hyperspectral processing patterns (memory-mapped arrays, chunked processing, interleave formats BIP/BIL/BSQ)
- `agent-architecture.md` - Target architecture with Data/Domain/Application/Presentation layers
- `agent-ui.md` - Tkinter UI patterns, keeping UI responsive with worker threads
- `agent-data.md` - File I/O patterns, lazy loading, metadata preservation
- `agent-feature.md` - Feature development workflow and quality checklist

**Key patterns:**
- Use vectorized NumPy operations for spectral processing
- Keep UI responsive - offload heavy computation
- Handle edge cases: empty ROI, single pixel, NaN values
- Preserve spectral fidelity - avoid lossy operations
- Use memory-mapped arrays for large datacubes
