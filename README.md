
# pyMelt mantle melting library

[![Python package](https://github.com/simonwmatthews/pyMelt/actions/workflows/python-package.yml/badge.svg)](https://github.com/simonwmatthews/pyMelt/actions/workflows/python-package.yml)
[![flake8](https://github.com/simonwmatthews/pyMelt/actions/workflows/flake8.yml/badge.svg)](https://github.com/simonwmatthews/pyMelt/actions/workflows/flake8.yml)
[![Documentation Status](https://readthedocs.org/projects/pymelt/badge/?version=latest)](https://pymelt.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/pyMelt.svg)](https://badge.fury.io/py/pyMelt)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/simonwmatthews/pyMelt/HEAD?labpath=docs%2Ftutorial%2Ftutorial1.ipynb)
[![DOI](https://zenodo.org/badge/259243892.svg)](https://zenodo.org/badge/latestdoi/259243892)

## Introduction

pyMelt is a python library for calculating the melting behaviour of mantle comprising multiple lithologies. The pyMelt library implements the melting equations developed by Phipps Morgan (2001), alongside many existing lherzolite and pyroxenite melting parameterisations.

Parameters that can be calculated:

* The geotherm for decompressing mantle
* Melt fractions for each lithology
* Crustal thickness for passive-upwelling at a mid-ocean ridge
* Crystallisation temperatures (following the method in Matthews et al., 2016)
* Magma flux at intra-plate settings
* Lava trace element concentrations

The `hydrousLithology` module allows hydrous melting to be approximated given any anhydrous lithology melting model using a modified version of the Katz et al. (2003) hydrous melting equations.

## Documentation
Full documentation, further information about the package, and a tutorial for getting started are
provided at [pymelt.readthedocs.io](http://pymelt.readthedocs.io).

## Installation
pyMelt is available on pip, and can be installed by running `pip install pyMelt` in a terminal.

## Run pyMelt on the cloud with myBinder
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/simonwmatthews/pyMelt/HEAD?labpath=docs%2Ftutorial%2Ftutorial1.ipynb)
You can use pyMelt and go through the tutorials right now without installing anything.

## pyMelt_MultiNest
pyMelt can be used in conjunction with the MultiNest algorithm (Feroz and Hobson, 2008; Feroz et al., 2009, 2013) via its python frontend, pyMultinest (Buchner et al., 2014). This permits the inversion of measured data (e.g. crystallisation temperature, crustal thickness) to obtain unknowns (e.g. potential temperature) via Bayesian inference. More details of the inversion methods are provided in Matthews et al., 2021.

For pyMelt_MultiNest to work, MultiNest and pyMultinest must be installed. The user is directed to the [pyMultinest installation instructions](https://johannesbuchner.github.io/PyMultiNest/) for further guidance.

Note that the pyMelt_MultiNest library is built for an old version of pyMelt and has not yet been updated.

## Citing pyMelt
If pyMelt enables or aids your research please cite the [pyMelt preprint](https://doi.org/10.31223/X5JP7X), hosted on EarthArXiv. To ensure reproducibility please also state the version of pyMelt that you used. The latest release is v1.960, which is archived in a Zenodo repository with DOI:
[![DOI](https://zenodo.org/badge/259243892.svg)](https://zenodo.org/badge/latestdoi/259243892)

A manuscript describing the module will be released as a preprint soon.

You should also cite the relevant publications for the pure-lithology melting models. If you use our models, you should cite:

Matthews, S., Wong, K., Shorttle, O., Edmonds, M., & Maclennan, J. (2021). Do olivine crystallization temperatures faithfully record mantle temperature variability?. Geochemistry, Geophysics, Geosystems, 22(4), e2020GC009157. https://doi.org/10.1029/2020GC009157
