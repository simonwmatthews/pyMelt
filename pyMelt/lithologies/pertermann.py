"""
================================
Pertermann and Hirschmann (2002)
================================

The pertermann module implements the G2 model.

"""

from pyMelt.lithology_classes import lithology as _lithology
from pyMelt.lithology_classes import default_properties as _default_properties

import numpy as _np


class g2(_lithology):
    """
    Implementation of the Pertermann and Hirschmann (2002) G2 melting model.
    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class, with
    values:

    - a:  Parameter used in calculating melt fraction.
    - b:  Parameter used in calculating melt fraction.
    - c:  Parameter used in liquidus definition.
    - d:  Parameter used in liquidus definition.
    - e:  Parameter used in solidus definition.
    - f:  Parameter used in solidus definition.

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    also be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: pyMelt.lithology_class.default_properties['CP']
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: pyMelt.lithology_class.default_properties['alphas']
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt.lithology_class.default_properties['alphaf']
        The thermal expansivity of the melt (1e-6 K-1)
    rhos :       float, default: pyMelt.lithology_class.default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt.lithology_class.default_properties['rhof']
        The density of the melt (kg m-3)
    DeltaS :     float, default: pyMelt.lithology_class.default_properties['DeltaS']
        The entropy of fusion J K-1 kg-1
    parameters : dict, default: parameters from Matthews et al. (2021)
        The model parameters described above
    """

    def __init__(self,
                 CP=_default_properties['CP'],
                 alphas=_default_properties['alphas'],
                 alphaf=_default_properties['alphaf'],
                 rhos=_default_properties['rhos'],
                 rhof=_default_properties['rhof'],
                 DeltaS=_default_properties['DeltaS'],
                 parameters={'a':     0.7368,
                             'b':     0.2632,
                             'c':   1175.000,
                             'd':    114.000,
                             'e':    920.000,
                             'f':    130.000
                             }
                 ):

        self.DeltaS = DeltaS
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.parameters = parameters

    def F(self, P, T, **kwargs):
        """
        Calculates melt fraction at a given pressure and temperature using:
            a*T'**2 + b*T',
        where T is the normalised temperature:
            (T-Tsolidus)/(T-Tliquidus).
        If P and T are below the the solidus, 0 is returned, if they are above the liquidus, 1 is
        returned.

        Parameters
        ----------
        P : float
            Pressure (GPa)
        T : float
            Temperature (degC)

        Returns
        -------
        float
            Melt fraction.
        """
        Tsol = self.TSolidus(P, **kwargs)
        Tliq = self.TLiquidus(P, **kwargs)
        if T < Tsol:
            F = 0.0
        elif T > Tliq:
            F = 1.0
        else:
            Tr = (T-Tsol)/(Tliq-Tsol)
            F = self.parameters['a']*Tr**2 + self.parameters['b']*Tr
        return F

    def TLiquidus(self, P, **kwargs):
        """
        Calculates the liquidus temperature, at a given pressure, using:
            c + d*P.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Liquidus temperature (degC).
        """
        Tliq = self.parameters['c'] + self.parameters['d']*P
        return Tliq

    def TSolidus(self, P, **kwargs):
        """
        Calculates the solidus temperature, at a given pressure, using:
            e + f*P.

        Parameters
        ----------
        P : float
            Pressure (GPa)

        Returns
        -------
        float
            Solidus temperature (degC).
        """
        Tsol = self.parameters['e'] + self.parameters['f']*P
        return Tsol

    def dTdF(self, P, T, **kwargs):
        """
        Calculates dT/dF(const. P) at a given pressure and temperature.

        Parameters
        ----------
        P : float
            Pressure (GPa).
        T : float
            Temperature (degC).

        Returns
        -------
        float
            dT/dF(const. P) (K)
        """
        Tsol = self.TSolidus(P, **kwargs)
        Tliq = self.TLiquidus(P, **kwargs)
        if T < Tsol:
            dTdF = _np.inf
        elif T > Tliq:
            dTdF = _np.inf
        else:
            dTdF = ((Tliq - Tsol) / (self.parameters['a'] * 2
                    * ((T - Tsol) / (Tliq - Tsol))+self.parameters['b']))

        return dTdF

    def dTdP(self, P, T, **kwargs):
        """
        Calculates dT/dP(const. F) at a given pressure and temperature.

        Parameters
        ----------
        P : float
            Pressure (GPa).
        T : float
            Temperature (degC).

        Returns
        -------
        float
            dTdP(const. F) (K GPa-1)
        """
        Tsol = self.TSolidus(P, **kwargs)
        Tliq = self.TLiquidus(P, **kwargs)
        dTdP = ((self.parameters['d'] - self.parameters['f'])
                * (T - Tsol)/(Tliq - Tsol) + self.parameters['f'])

        return dTdP
