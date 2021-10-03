"""
=========
Chemistry
=========

The chemistry module provides the base classes for defining chemical elements/species for
inclusion in pyMelt calculations, alongside default implementations.
"""

import numpy as _np


class species(object):
    """
    The species base class. The class contains all the parameters and methods required to
    calculate the concentration of chemical species in melts, given T, P, F, etc.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.

    """

    def __init__(self, name, c0, **kwargs):
        self.name = name
        self.c0 = c0

    def instantaneous_melt(self, state):
        """
        Returns the concentration of the species in the instantaneous melt for the specified state.
        This function should be redefined according to the chemical model being used. If the model
        does not predict instantaneous melts it need not be redefined.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        Returns
        -------
        float
            The concentration of the species in the melt.
        """
        return _np.nan

    def accumulated_melt(self, state):
        """
        Returns the concentration of the species in the accumulated melt at the specified state.
        This function should be redefined according to the chemical model being used. If the model
        does not predict accumulated melts it need not be redefined.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        Returns
        -------
        float
            The concentration of the species in the melt.
        """
        return _np.nan


class BatchSpecies(species):
    """
    Implementation of the species class for batch melting with a constant partition coefficient.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.
    D : float
        The partition coefficient
    """

    def __init__(self, name, c0, D, **kwargs):
        self.name = name
        self.c0 = c0
        self.D = D

    def accumulated_melt(self, state):
        """
        Returns the concentration in the melt during batch melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """
        return self.c0 / (self.D * (1 - state['F']) + state['F'])


class ContinuousSpecies(species):
    """
    Implementation of the species class for batch melting with a constant partition coefficient.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.
    D : float
        The partition coefficient
    """

    def __init__(self, name, c0, D, phi=0.005, **kwargs):
        self.name = name
        self.c0 = c0
        self.D = D
        self.phi = phi

    def instantaneous_melt(self, state):
        """
        Returns the instantaneous concentration in the melt during near-fractional (continuous)
        melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """

        exponent = (1 - self.phi) * (1 - self.D) / ((1 - self.phi) * self.D + self.phi)

        return self.c0 / ((1 - self.phi) * self.D + self.phi) * (1 - state.F)**exponent

    def accumulated_melt(self, state):
        """
        Returns the concentration in the melt during batch melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """

        Deff = (1 - self.phi) * self.D + self.phi

        return self.c0 / state.F * (1 - (1 - state.F)**(1 / Deff))


workman05_ddm = {'Rb': 0.05,
                 'Ba': 0.563,
                 'Th': 0.0079,
                 'U': 0.0032,
                 'Nb': 0.1485,
                 'Ta': 0.0096,
                 'La': 0.192,
                 'Ce': 0.550,
                 'Pb': 0.018,
                 'Pr': 0.107,
                 'Nd': 0.581,
                 'Sr': 7.664,
                 'Zr': 5.082,
                 'Hf': 0.157,
                 'Sm': 0.239,
                 'Eu': 0.096,
                 'Ti': 716.3,
                 'Gd': 0.358,
                 'Tb': 0.070,
                 'Dy': 0.505,
                 'Ho': 0.115,
                 'Y': 3.328,
                 'Er': 0.348,
                 'Yb': 0.365,
                 'Lu': 0.058
                 }

workman05_D = {'Rb': 1e-5,
               'Ba': 0.00012,
               'Th': 0.001,
               'U': 0.0011,
               'Nb': 0.0034,
               'Ta': 0.0034,
               'La': 0.01,
               'Ce': 0.022,
               'Pb': 0.014,
               'Pr': 0.027,
               'Nd': 0.031,
               'Sr': 0.025,
               'Zr': 0.033,
               'Hf': 0.035,
               'Sm': 0.045,
               'Eu': 0.050,
               'Ti': 0.058,
               'Gd': 0.056,
               'Tb': 0.068,
               'Dy': 0.079,
               'Ho': 0.084,
               'Y': 0.088,
               'Er': 0.097,
               'Yb': 0.115,
               'Lu': 0.120
               }

stracke03_bsic = {'Rb': 0.57,
                  'Ba': 6.59,
                  'Th': 0.088,
                  'U': 0.027,
                  'Nb': 1.95,
                  'Ta': 0.124,
                  'La': 1.68,
                  'Ce': 5.89,
                  'Pb': 0.09,
                  'Pr': 0.0,
                  'Nd': 7.45,
                  'Sr': 81.0,
                  'Zr': 64.0,
                  'Hf': 1.78,
                  'Sm': 2.69,
                  'Eu': 1.04,
                  'Ti': 7735.0,
                  'Gd': 4.03,
                  'Tb': 0.0,
                  'Dy': 5.01,
                  'Ho': 0.0,
                  'Y': 28.5,
                  'Er': 3.13,
                  'Yb': 2.99,
                  'Lu': 0.45}
