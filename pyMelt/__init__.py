"""
======
pyMelt
======

A python package for performing mantle melting calculations.

Developed and maintained by Simon Matthews (simonm@hi.is) and Kev Wong.

This module facilitates calculations for melting of heterogeneous mantle, based on empirical
parameterisations of the melting functions. The calculations assume complete thermal equilibrium
and complete chemical disequilibrium between adjacent mantle lithologies. Various tools are
provided for performing decompression melting calculations.
"""

__version__ = "0.5.1"
__author__ = "Simon Matthews and Kevin Wong"

import pyMelt.lithologies
# import pyMelt.meltingcolumn_classes
import pyMelt.mantle_class
from pyMelt.lithology_classes import hydrous_lithology


class Mantle(pyMelt.mantle_class.mantle):
    __doc__ = pyMelt.mantle_class.mantle.__doc__
