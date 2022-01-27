"""
==================
McKenzie and Bickle (1988)
==================

Implementation of the anhydrous melting model presented by McKenzie and Bickle (1988).

"""

from pyMelt.lithology_classes import lithology as _lithology
from pyMelt.lithology_classes import default_properties as _default_properties
import numpy as np
from scipy.optimize import fsolve
from scipy.special import expit


class lherzolite(_lithology):
    """
    Implementation of the McKenzie and Bickle (1988) garnet peridotite melting model.
    As this parameterisation provides pressure as a function of solidus temperature,
    scipt.optimize.fsolve and scipy.special.expit are required to numerically find
    solidus temperature as a function of pressure.
    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class, with
    values:
    - A1:  Parameter used in solidus definition.
    - A2:  Parameter used in solidus definition.
    - A3:  Parameter used in solidus definition.
    - A4:  Parameter used in solidus definition.
    - B1:  Parameter used in liquidus definition.
    - B2:  Parameter used in liquidus definition.
    - B3:  Parameter used in liquidus definition.
    - B4:  Parameter used in liquidus definition.
    - a0:  Parameter used in calculating melt fraction.
    - a1:  Parameter used in calculating melt fraction.
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
    parameters : dict, default: parameters from McKenzie and Bickle (1988)
        The model parameters described above
    """

    def __init__(self,
                 CP=_default_properties['CP'],
                 alphas=_default_properties['alphas'],
                 alphaf=_default_properties['alphaf'],
                 rhos=_default_properties['rhos'],
                 rhof=_default_properties['rhof'],
                 DeltaS=_default_properties['DeltaS'],
                 parameters={'A1': 1100.0,
                             'A2': 136.0,
                             'A3': 4.968e-4,
                             'A4': 1.2e-2,
                             'B1': 1736.0,
                             'B2': 4.343,
                             'B3': 180.0,
                             'B4': 2.2169,
                             'a0': 0.4256,
                             'a1': 2.988
                             }
                 ):
        self.DeltaS = DeltaS
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.parameters = parameters

    def TSolidus(self, P):
        """
        Calculates the solidus temperature, at a given pressure, using Equation 18 of
        McKenzie and Bickle (1988). Requires scipy.optimize.fsolve to solve for pressure.
        As fsolve can only solve for a single pressure TSolidus first detects whether P is
        given as a list, tuple, numpy.ndarray (e.g. numpy.linspace), or scalar prior to
        calculating the solidus pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa)
        Returns
        -------
        float
            Solidus temperature (degC).
        """
        TSolidus_initial_guess = 1300.0

        def func(TSolidus):
            _func = ((TSolidus - self.parameters['A1']) / self.parameters['A2']
                     + self.parameters['A3'] * np.exp(self.parameters['A4']
                     * (TSolidus - self.parameters['A1']))) - Pressure
            return _func

        if isinstance(P, (list, tuple, np.ndarray)) is True:
            TSolidus_solution = []
            for i in range(len(P)):
                Pressure = P[i]
                TSolidus_solution.append(fsolve(func, TSolidus_initial_guess, maxfev=50000)[0])
        else:
            Pressure = P
            TSolidus_solution = fsolve(func, TSolidus_initial_guess, maxfev=50000)[0]
        return TSolidus_solution

    def TLiquidus(self, P):
        """
        Calculates the liquidus temperature, at a given pressure, using Equation 19 of
        McKenzie and Bickle (1988).

        Parameters
        ----------
        P : float
            Pressure (GPa)
        Returns
        -------
        float
            Solidus temperature (degC).
        """
        TLiquidus = (self.parameters['B1'] + self.parameters['B2'] * P + self.parameters['B3']
                     * np.arctan(P / self.parameters['B4']))
        return TLiquidus

    def _dTdPSolidus(self, P):
        """
        Returns the solidus temperature gradient at any given pressure. Requires
        scipy.optimize.fsolve to solve for pressure, and scipy.special.expit to
        avoid exponential overflows that can occur with numpy.exp.

        Parameters
        ----------
        P : float
            Pressure (GPa).
        Returns
        -------
        float
            Solidus temperaure gradient (degC/GPa)
        """
        TSolidus = self.TSolidus(P)
        dTdPSolidus = (1 / (1 / self.parameters['A2'] + self.parameters['A3']
                       * self.parameters['A4'] * expit(self.parameters['A4']
                       * (TSolidus - self.parameters['A1']))))
        return dTdPSolidus

    def _dTdPLiquidus(self, P):
        """
        Returns the liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).
        Returns
        -------
        float
            Solidus temperaure gradient (degC/GPa)
        """
        dTdPLiquidus = (self.parameters['B2'] + self.parameters['B3']
                        / (1 + (P / self.parameters['B4'])**2))
        return dTdPLiquidus

    def _RescaledT(self, T, P):
        """
        Calculates the rescaled temperature defined by Equation 20 of
        McKenzie and Bickle (1988).

        Parameters
        ----------
        T : float
            Temperature (degC).
        P : float
            Pressure (GPa).
        Returns
        -------
        float
            Rescaled Temperature (dimensionless).
        """
        TSolidus = self.TSolidus(P)
        TLiquidus = self.TLiquidus(P)
        RescaledT = (T - (TSolidus + TLiquidus) / 2) / (TLiquidus - TSolidus)
        return RescaledT

    def F(self, P, T):
        """
        Wrapper for the melt fraction functions, defined by Equations 21 and 22
        of McKenzie and Bickle (1988). If T and P are below the solidus, returns 0,
        if they are above the liquidus, returns 1. Otherwise determines F.

        Parameters
        ----------
        P : float
            Pressure (GPa).
        T : float
            Temperature (degC).
        Returns
        -------
        float
            Melt fraction.
        """
        TSolidus = self.TSolidus(P)
        TLiquidus = self.TLiquidus(P)
        RescaledT = self._RescaledT(T, P)
        if T > TLiquidus:
            F = 1.0
        elif T < TSolidus:
            F = 0.0
        else:
            F = (RescaledT + (RescaledT**2 - 0.25) * (self.parameters['a0']
                 + self.parameters['a1'] * RescaledT) + 0.5)
        return F

    def dTdF(self, P, T):
        """
        Calculates dT/dF(const. P). First calculates the melt fraction. If F is
        zero, returns _np.inf. If F is 1, returns _np.inf. Otherwise uses the
        appropriate expression for melt fraction.

        Parameters
        ----------
        P : float
            Pressure (GPa)
        T : float
            Temperature (degC)
        Returns
        -------
        float
            dT/dF(const. P) (K).
        """
        TSolidus = self.TSolidus(P)
        TLiquidus = self.TLiquidus(P)
        RescaledT = self._RescaledT(T, P)
        if T < TSolidus:
            dTdF = np.inf
        elif T > TLiquidus:
            dTdF = np.inf
        else:
            dFdT = (((1 - 0.25 * self.parameters['a1']) + 3 * self.parameters['a0'] * RescaledT
                     + 3 * self.parameters['a1'] * RescaledT**2) / (TLiquidus - TSolidus))
            dTdF = 1 / dFdT
        return dTdF

    def dTdP(self, P, T):
        """
        Calculates dT/dP(const. F). First calculates F, then chooses the
        appropriate expression for melt fraction.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).
        Returns
        -------
        float
            dT/dP(const. F) (K GPa-1).
        """
        dTdPSolidus = self._dTdPSolidus(P)
        dTdPLiquidus = self._dTdPLiquidus(P)
        RescaledT = self._RescaledT(T, P)
        dTdP = RescaledT * (dTdPLiquidus - dTdPSolidus) + 0.5 * (dTdPSolidus + dTdPLiquidus)
        return dTdP
