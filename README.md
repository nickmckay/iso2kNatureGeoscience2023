[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8327340.svg)](https://doi.org/10.5281/zenodo.8327340)


# Code to reproduce key analyses and figures from Konecky et al,. (2023) "Globally coherent water cycle response to temperature change during the past two millennia"

This repository includes the R, Matlab and Julia code used to calculate the results shown in the figures. For some analyses/figures the plotting code is also included. 

## Figures

The codesets for each figure are included as subdirectories of the "Figures" directory. 

## Data

Where possible, we included code to download the necessary data directly in the scripts. These datasets are also available in the "data" directory. The major exception to this are the CESM data, which are too large to store here and are available on the NCAR campaign storage facility.

## Outputs

Where the included codesets create new data products, they are included in the "data/outputs" directory.

## Utilities

The "utilities" directory includes supplemental codebases referenced in the figure codesets that users may need to reference.

### Spatial Data

Data and metadata needed to recreate some of the maps are included in the "spatialData" subdirectory

### Markov Block Bootstrap

Install using the Julia [package manager](https://pkgdocs.julialang.org/v1/): 

```julia
import Pkg
Pkg.add(path="https://github.com/nickmckay/iso2kNatureGeoscience2023", subdir="utilities/Juliapkgs/MarkovBlocks.jl")
```

## How to cite this repository

McKay, Nicholas; Konecky, Bronwen; Falster, Georgina; Stevenson, Samantha and  Fischer, Matt. Code to reproduce key analyses and figures from Konecky et al,. (2023) "Globally coherent water cycle response to temperature change during the past two millennia". 2023. Zenodo, doi:10.5281/zenodo.8327340.
