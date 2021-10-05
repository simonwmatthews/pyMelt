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
                  'Pr': _np.nan,
                  'Nd': 7.45,
                  'Sr': 81.0,
                  'Zr': 64.0,
                  'Hf': 1.78,
                  'Sm': 2.69,
                  'Eu': 1.04,
                  'Ti': 7735.0,
                  'Gd': 4.03,
                  'Tb': _np.nan,
                  'Dy': 5.01,
                  'Ho': _np.nan,
                  'Y': 28.5,
                  'Er': 3.13,
                  'Yb': 2.99,
                  'Lu': 0.45}

olv_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0005,
         'Ce': 0.0005,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.00042,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.0011,
         'Eu': 0.0016,
         'Ti': _np.nan,
         'Gd': 0.0011,
         'Tb': 0.0015,
         'Dy': 0.0027,
         'Ho': 0.0016,
         'Y': _np.nan,
         'Er': 0.013,
         'Yb': 0.020,
         'Lu': 0.020,
         }

opx_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0031,
         'Ce': 0.0040,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.01200,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.0200,
         'Eu': 0.0130,
         'Ti': _np.nan,
         'Gd': 0.0130,
         'Tb': 0.0190,
         'Dy': 0.0110,
         'Ho': 0.0065,
         'Y': _np.nan,
         'Er': 0.045,
         'Yb': 0.080,
         'Lu': 0.120,
         }

cpx_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0490,
         'Ce': 0.0800,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.17800,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.2930,
         'Eu': 0.3350,
         'Ti': _np.nan,
         'Gd': 0.3500,
         'Tb': 0.4030,
         'Dy': 0.4000,
         'Ho': 0.4270,
         'Y': _np.nan,
         'Er': 0.420,
         'Yb': 0.400,
         'Lu': 0.376,
         }

grt_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0010,
         'Ce': 0.0050,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.05200,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.2500,
         'Eu': 0.4960,
         'Ti': _np.nan,
         'Gd': 0.84800,
         'Tb': 1.4770,
         'Dy': 2.2000,
         'Ho': 3.3150,
         'Y': _np.nan,
         'Er': 4.400,
         'Yb': 6.600,
         'Lu': 7.100,
         }

spn_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0100,
         'Ce': 0.0100,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.0100,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.0100,
         'Eu': 0.0100,
         'Ti': _np.nan,
         'Gd': 0.0100,
         'Tb': 0.0100,
         'Dy': 0.0100,
         'Ho': 0.0100,
         'Y': _np.nan,
         'Er': 0.0100,
         'Yb': 0.0100,
         'Lu': 0.0100,
         }

plg_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.2700,
         'Ce': 0.200,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.1400,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.1100,
         'Eu': 0.7300,
         'Ti': _np.nan,
         'Gd': 0.0660,
         'Tb': 0.0600,
         'Dy': 0.0550,
         'Ho': 0.0480,
         'Y': _np.nan,
         'Er': 0.0100,
         'Yb': 0.031,
         'Lu': 0.0250,
         }
