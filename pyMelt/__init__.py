"""pyMelt - a python package for performing mantle melting calculations.

Developed and maintained by Simon Matthews (simonm@hi.is) and Kev Wong.

This module facilitates calculations for melting of heterogeneous mantle, based on empirical
parameterisations of the melting functions. The calculations assume complete thermal equilibrium
and complete chemical disequilibrium between adjacent mantle lithologies. Various tools are
provided for performing decompression melting calculations.
"""

import pyMelt.lithology_class
import pyMelt.lithologies
import pyMelt.mantle_class


class Mantle(pyMelt.mantle_class.mantle):
    """
    The mantle class consists of one or more lithology classes, in a particular proportion. The
    mantle class contains the methods used for doing melting calculations and for calculating the
    properties of a heterogeneous mantle.

    Parameters
    ----------
    lithologies :    list of lithology objects
        A list of lithology instances.
    proportions :    list or numpy.array of floats
        The mass ratios of the lithologies, doesn't need to be normalised.
    names :  list of strings or None, default: None
        The names of the lithologies. If False, default names will be chosen.

    Attributes
    ----------
    number_lithologies : int
        the number of lithologies in the mantle class
    CP :     list of floats
        the heat capacities of the lithologies
    alphaf : list of floats
        the thermal expansion coefficients of the melts produced by each lithology.
    alphas : list of floats
        the thermal expansion coefficients of each lithology
    rhof :   list of floats
        the densities of the melts produced by each lithology
    rhos :   list of floats
        the densities of each lithology
    DeltaS : list of floats
        the entropy change on melting of each lithology.

    """
    pass


class Lithology(pyMelt.lithology_class.lithology):
    """
    Lithology base class. This class contains all the parameters and methods required to calculate
    the melting behaviour of a single mantle component.

    Attributes
    ----------
    DeltaS:     float
        Entropy of fusion. (J kg-1 K-1). Default is 300.0.
    CP:     float
        Heat capacity (J Kg-1 K-1). Default is 1000.0.
    alphas:     float
        Thermal expansivity of the solid (K-1). Default is 40.0.
    alphaf:     float
        Thermal expansivity of the melt (K-1). Default is 68.0.
    rhos:   float
        Density of the solid (g cm-3). Default is 3.3.
    rhof:   float
        Density of the melt (g cm-3). Default is 2.9.

    """
    pass
