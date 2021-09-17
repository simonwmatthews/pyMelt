
# pyMelt mantle melting library
S. Matthews (University of Iceland) and K. Wong (University of Leeds)

## Introduction

pyMelt is a python package for calculating the melting behaviour of mantle comprising multiple lithologies. The module implements the melting equations developed by Phipps Morgan (2001) to calculate the melting behaviour of mantle comprising any chosen lithology.

Currently supported calculations:

* Adiabatic decompression melting
* Isobaric melting

Parameters that can be calculated:

* The geotherm for decompressing mantle
* Melt fractions for each lithology
* Crustal thickness for passive-upwelling at a mid-ocean ridge
* Crystallisation temperatures (following the method in Matthews et al., 2016)

## Documentation
Full documentation, further information about the package, and a tutorial for getting started are
provided at [pymelt.readthedocs.io](http://pymelt.readthedocs.io).

## Installation
pyMelt is available on pip, and can be installed by running `pip install pyMelt` in a terminal.

## pyMelt_MultiNest
pyMelt can be used in conjunction with the MultiNest algorithm (Feroz and Hobson, 2008; Feroz et al., 2009, 2013) via its python frontend, pyMultinest (Buchner et al., 2014). This permits the inversion of measured data (e.g. crystallisation temperature, crustal thickness) to obtain unknowns (e.g. potential temperature) via Bayesian inference. More details of the inversion methods are provided in Matthews et al., 2021.

For pyMelt_MultiNest to work, MultiNest and pyMultinest must be installed. The user is directed to the [pyMultinest installation instructions](https://johannesbuchner.github.io/PyMultiNest/) for further guidance.

## Citing pyMelt
If pyMelt enables or aids your research please cite the release you used. The latest release is v1.915 and has the doi:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5513675.svg)](https://doi.org/10.5281/zenodo.5513675)

You should also cite the relevant publications for the pure-lithology melting models. If you use our models, you should cite:

Matthews, S., Wong, K., Shorttle, O., Edmonds, M., & Maclennan, J. (2021). Do olivine crystallization temperatures faithfully record mantle temperature variability?. Geochemistry, Geophysics, Geosystems, 22(4), e2020GC009157. https://doi.org/10.1029/2020GC009157
