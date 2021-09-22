"""
======================
Shorttle et al. (2014)
======================

Implementation of the new lithologies in Shorttle et al. (2014).

"""

from pyMelt.lithology_classes import lithology as _lithology
from pyMelt.lithology_classes import default_properties as _default_properties

import numpy as _np


class kg1(_lithology):
    """
    Implementation of the KG1 parameterisation by Shorttle et al. (2014).

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class, with
    values:

    - A1:     Constant used in solidus expression.
    - A2:     Constant used in solidus expression.
    - A3:     Constant used in solidus expression.
    - B1:     Constant used in cpx-out expression.
    - B2:     Constant used in cpx-out expression.
    - B3:     Constant used in cpx-out expression.
    - C1:     Constant used in liquidus expression.
    - C2:     Constant used in liquidus expression.
    - C3:     Constant used in liquidus expression.
    - a:      Constant used in cpx-present melt fraction expression.
    - b:      Constant used in cpx-present melt fraction expression.
    - c:      Constant used in cpx-absent melt fraction expression.
    - d:      Constant used in cpx-absent melt fraction expression.
    - alpha:  Exponent used in the cpx-present melt fraction expression.
    - beta:   Exponent used in the cpx-absent melt fraction expression.

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
    parameters : dict, default: parameters from Shorttle et al. (2014)
        The model parameters described above
    """
    def __init__(self,
                 CP=_default_properties['CP'],
                 alphas=_default_properties['alphas'],
                 alphaf=_default_properties['alphaf'],
                 rhos=_default_properties['rhos'],
                 rhof=_default_properties['rhof'],
                 DeltaS=_default_properties['DeltaS'],
                 parameters={'A1':    1095.4,
                             'A2':     124.1,
                             'A3': - 4.7,
                             'B1':    1179.6,
                             'B2':     157.2,
                             'B3': - 11.1,
                             'C1':    1780.0,
                             'C2':      45.0,
                             'C3': - 2.0,
                             'a':        0.3187864,
                             'b':        0.4154,
                             'c':        0.7341864,
                             'd':        0.2658136,
                             'alpha':    2.0,
                             'beta':     1.5
                             }
                 ):

        self.DeltaS = DeltaS
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.parameters = parameters

    def TSolidus(self, P, **kwargs):
        """
        Returns solidus temperature at a given pressure. T = A1 + A2*P + A3*P**2.

        Parameters
        ----------
        P : float, or list like
            Pressure in GPa.

        Returns
        -------
        float, or list like
            Solidus Temperature in degC.
        """
        T = self.parameters['A1'] + self.parameters['A2']*P + self.parameters['A3']*P**2
        return T

    def TLiquidus(self, P, **kwargs):
        """
        Returns liquidus temperature at a given pressure. T = C1 + C2*P + C3*P**2.

        Parameters
        ----------
        P :  float, or list like
            Pressure in GPa.

        Returns
        -------
        float, or list like
            Liquidus temperature in degC.
        """
        T = self.parameters['C1'] + self.parameters['C2']*P + self.parameters['C3']*P**2
        return T

    def dTdP(self, P, T, **kwargs):
        """
        Returns dT/dP (constant F) at a given pressure and temperature.

        Parameters
        ----------
        P : float
            Pressure in GPa.
        T : float
            Temperature in degC.

        Returns
        -------
        float
            dT/dP (constant F) in K GPa-1.
        """
        Tsol = self.TSolidus(P, **kwargs)
        Tliq = self.TLiquidus(P, **kwargs)
        Tcpx = self._TCpxOut(P, **kwargs)

        if T < Tcpx:
            dTdP = (-(-(T - Tsol) / (Tcpx - Tsol) *
                    (self.parameters['B2'] + 2*self.parameters['B3']*P -
                     self.parameters['A2'] - self.parameters['A3']*2*P) -
                    self.parameters['A2'] - self.parameters['A3']*2*P))
        else:
            dTdP = (-(-(T - Tcpx) / (Tliq - Tcpx) *
                    (self.parameters['C2'] + self.parameters['C3']*2*P -
                     self.parameters['B2'] - self.parameters['B3']*2*P) -
                    self.parameters['B2'] - 2*self.parameters['B3']*P))

        return dTdP

    def dTdF(self, P, T, **kwargs):
        """
        Returns dT/dF (constant P) at a given pressure and temperature. If below
        the solidus, or above the liquidus, _np.inf is returned.

        Parameters
        ----------
        P : float
            Pressure in GPa.
        T : float
            Temperature in degC.

        Returns
        -------
        float
            dT/dF (constant P) in K.
        """
        Tsol = self.TSolidus(P, **kwargs)
        Tliq = self.TLiquidus(P, **kwargs)
        Tcpx = self._TCpxOut(P, **kwargs)

        if T < Tsol:
            dTdF = _np.inf
        elif T < Tcpx:
            dTdF = ((Tcpx-Tsol) /
                    (self.parameters['b'] + self.parameters['a'] *
                     self.parameters['alpha'] *
                     ((T-Tsol)/(Tcpx-Tsol))**(self.parameters['alpha']-1)))
        elif T < Tliq:
            dTdF = ((Tliq - Tcpx)/(self.parameters['d'] * self.parameters['beta'] *
                    ((T-Tcpx)/(Tliq-Tcpx))**(self.parameters['beta']-1)))
        else:
            dTdF = _np.inf

        return dTdF

    def F(self, P, T, **kwargs):
        """
        Returns melt fraction at a given pressure and temperature. If below the
        solidus, returns 0. If above the liquidus, returns 1.

        Prior to cpx exhaustion:
        F = a*(Tr)**alpha _ b*Tr
        Tr = (T-Tsol)/(Tliq-Tsol)

        After cpx exhaustion:
        F = d*(Tr)**beta + c
        Tr = (T-Tcpx)/(Tliq-Tcpx)

        Parameters
        ----------
        P : float
            Pressure in GPa.
        T : float
            Temperature in degC.

        Returns
        -------
        float
            Melt fraction between 0 and 1.
        """
        Tsol = self.TSolidus(P, **kwargs)
        Tliq = self.TLiquidus(P, **kwargs)
        Tcpx = self._TCpxOut(P, **kwargs)

        if T < Tsol:
            F = 0.0
        elif T > Tliq:
            F = 1.0
        elif T < Tcpx:
            Tf = (T-Tsol)/(Tcpx-Tsol)
            F = self.parameters['a'] * Tf**self.parameters['alpha'] + self.parameters['b']*Tf
        else:
            Tf = (T-Tcpx)/(Tliq-Tcpx)
            F = self.parameters['d'] * Tf**self.parameters['beta'] + self.parameters['c']
        return F

    def _TCpxOut(self, P, **kwargs):
        """
        Returns the temperature of cpx-exhaustion at a given pressure. T = B1 + B2*P + B3*P**2.

        Parameters
        ----------
        P :float, or list like
            Pressure in GPa.

        Returns
        -------
        float, or list like
            Cpx-exhaustion temperature in degC.
        """
        T = self.parameters['B1'] + self.parameters['B2']*P + self.parameters['B3']*P**2
        return T


class harzburgite(_lithology):
    """
    Material that does not melt, i.e. Harzburgite in Shorttle et al. (2014) and
    Matthews et al. (2016). Default thermodynamic constants are those used by
    Katz et al. (2003).

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: pyMelt._default_properties['CP']
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: pyMelt._default_properties['alphas']
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt._default_properties['alphaf']
        Melt thermal expansivity, not used, here for consistency.
    rhos :       float, default: pyMelt._default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt._default_properties['rhof']
        Melt density, not used, here for consistency.
    DeltaS :     float, default: pyMelt._default_properties['DeltaS']
        The entropy of fusion, not used, here for consistency.
    parameters : dict, default: {}
        This model does not use any parameters, here for consistency.
    """
    def __init__(self,
                 CP=_default_properties['CP'],
                 alphas=_default_properties['alphas'],
                 alphaf=_default_properties['alphaf'],
                 rhos=_default_properties['rhos'],
                 rhof=_default_properties['rhof'],
                 DeltaS=_default_properties['DeltaS'],
                 parameters={}
                 ):
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.DeltaS = DeltaS
        self.parameters = parameters

    def F(self, P, T, **kwargs):
        """
        Melt Fraction. Returns 0.0.

        Parameters
        ----------
        P : float
            Pressure. There to maintain consistancy within lithology definitions.
        T : float
            Temperature. There to maintain consistancy within lithology definitions.

        Returns
        -------
        float
            Melt fraction will always be 0.0.
        """
        return 0.0

    def dTdF(self, P, T, **kwargs):
        """
        dTdF(constP). Returns _np.inf.

        Parameters
        ----------
        P : float
            Pressure. There to maintain consistancy within lithology definitions.
        T : float
            Temperature. There to maintain consistancy within lithology definitions.

        Returns
        -------
        numpy.inf
            The value will always be infinite.
        """
        return _np.inf

    def dTdP(self, P, T, **kwargs):
        """
        dTdP(constF). Returns 0.0.

        Parameters
        ----------
        P : float
            Pressure. There to maintain consistancy within lithology definitions.
        T : float
            Temperature. There to maintain consistancy within lithology definitions.

        Returns
        -------
        float
            The value will always be 0.0.
        """
        return 0.0

    def TSolidus(self, P, **kwargs):
        """
        Solidus temperature. Returns _np.inf.

        Parameters
        ----------
        P : float
            Pressure. There to maintain consistancy within lithology definitions.

        Returns
        -------
        numpy.inf
            The value will always be infinite.
        """
        if isinstance(P, list) or isinstance(P, _np.ndarray):
            return _np.array([_np.inf]*len(P))
        else:
            return _np.inf

    def TLiquidus(self, P, **kwargs):
        """
        Liquidus temperature. Returns _np.inf

        Parameters
        ----------
        P : float
            Pressure. There to maintain consistancy within lithology definitions.

        Returns
        -------
        numpy.inf
            The value will always be infinite.
        """
        if isinstance(P, list) or isinstance(P, _np.ndarray):
            return _np.array([_np.inf]*len(P))
        else:
            return _np.inf
