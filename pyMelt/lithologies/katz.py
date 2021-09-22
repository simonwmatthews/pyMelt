"""
==================
Katz et al. (2003)
==================

Implementation of the anhydrous melting model presented by Katz et al. (2003).

"""

from pyMelt.lithology_classes import lithology as _lithology
from pyMelt.lithology_classes import default_properties as _default_properties

import numpy as _np


class lherzolite(_lithology):
    """
    Implementation of the Katz et al. (2003) anhydrous lherzolite melting model.

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class.

    - Mcpx:   Mass fraction of cpx in the source. Controls the transition to low-productivity
      harzburgite-type melting.
    - A1:     Parameter used to define solidus.
    - A2:     Parameter used to define solidus.
    - A3:     Parameter used to define solidus.
    - B1:     Parameter used to define lherzolite-liquidus.
    - B2:     Parameter used to define lherzolite-liquidus.
    - B3:     Parameter used to define lherzolite-liquidus.
    - C1:     Parameter used to define liquidus.
    - C2:     Parameter used to define liquidus.
    - C3:     Parameter used to define liquidus.
    - beta1:  Parameter used to calculate melt fraction during cpx-present melting.
    - beta2:  Parameter used to calculate melt fraction during cpx-absent melting.
    - r1:     Parameter used to define cpx reaction coefficient.
    - r2:     Parameter used to define cpx reaction coefficient.

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
                 parameters={'Mcpx':  0.15,
                             'A1': 1085.70,
                             'A2':  132.9,
                             'A3': - 5.1,
                             'B1': 1475.0,
                             'B2':   80.0,
                             'B3': - 3.2,
                             'C1': 1780.0,
                             'C2':   45.0,
                             'C3': - 2.0,
                             'beta1': 1.5,
                             'beta2': 1.5,
                             'r1':    0.5,
                             'r2':    0.08
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
        Returns the temperature of the solidus at any given pressure. Eqn(4).

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Solidus temperature (degC).
        """
        TSolidus = (self.parameters['A1'] + self.parameters['A2']*P + self.parameters['A3']*(P**2))
        return TSolidus

    def TLiquidus(self, P, **kwargs):
        """
        Returns the temperature of the liquidus at any given pressure. Eqn(10).

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Liquidus temperature (degC).
        """
        TLiquidus = self.parameters['C1'] + self.parameters['C2']*P + self.parameters['C3']*(P**2)
        return TLiquidus

    def F(self, P, T, **kwargs):
        """
        Wrapper for the melt fraction functions. If T and P are below the solidus, returns 0, if
        they are above the liquidus, returns 1. If below the temperature of cpx-exhaustion, calls
        the Fcpx function, otherwise calls the Fopx function.

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
        if T > self.TLiquidus(P, **kwargs):
            F = 1.0
        elif T < self.TSolidus(P, **kwargs):
            F = 0.0
        elif T < self._TcpxOut(P, **kwargs):
            F = self._Fcpx(T, P, **kwargs)
        else:
            F = self._Fopx(T, P, **kwargs)

        return F

    def dTdF(self, P, T, **kwargs):
        """
        Calculates dT/dF(const. P). First calculates the melt fraction. If F is zero, returns
        _np.inf. If F is 1, returns _np.inf. Otherwise uses the appropriate expressions for cpx
        present or absent melting.

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
        F = self.F(P, T, **kwargs)
        if F == 0:
            # If no melt fraction the derivative is zero. Prevents division by zero.
            dTdF = _np.inf
        elif F < self._FcpxOut(P, **kwargs):
            dTdF = (((1/self.parameters['beta1']))
                    * (self._TLherzLiquidus(P, **kwargs) - self.TSolidus(P, **kwargs))
                    * (F**((1-self.parameters['beta1'])/self.parameters['beta1'])))
        elif F < 1.0:
            dTdF = (((1/self.parameters['beta2']))
                    * (self.TLiquidus(P, **kwargs) - self._TcpxOut(P, **kwargs))
                    * (F**((1-self.parameters['beta2'])/self.parameters['beta2'])))
        else:
            dTdF = _np.inf
        return dTdF

    def dTdP(self, P, T, **kwargs):
        """
        Calculates dT/dP(const. F). First calculates F, then chooses the
        appropriate expression for cpx present or absent melting.

        Parameters
        ----------
        P : float
            Pressure (GPa).
        T : float
            Temperature (degC).

        Returns
        -------
        float
            dT/dP(const. F) (K GPa-1).
        """

        F = self.F(P, T, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        dTdPSolidus = self._dTdPSolidus(P, **kwargs)
        TLiquidus = self.TLiquidus(P, **kwargs)
        dTdPLiquidus = self._dTdPLiquidus(P, **kwargs)
        dTdPLherzLiquidus = self._dTdPLherzLiquidus(P, **kwargs)
        TcpxOut = self._TcpxOut(P, **kwargs)
        dTdPcpxOut = self._dTdPcpxOut(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        dFdPcpxOut = self._dFdPcpxOut(P, **kwargs)

        if F == 0:
            dTdP = self.alphas/self.rhos/self.CP
        elif F < self.FcpxOut(P, **kwargs):
            dTdP = (((F**(1/self.parameters['beta1']))
                    * (dTdPLherzLiquidus-dTdPSolidus)) + dTdPSolidus)
        elif F < 1.0:
            Trel = (T - TcpxOut)/(TLiquidus - TcpxOut)
            dTdP = ((TLiquidus - TcpxOut)/(1 - FcpxOut)
                    * (1/self.parameters['beta2'])*Trel**(1-self.paraeters['beta2'])
                    * dFdPcpxOut * (Trel**self.parameters['beta2'] - 1)
                    + dTdPcpxOut + Trel*(dTdPLiquidus - dTdPcpxOut))
        else:
            dTdP = self.alphaf/self.rhof/self.CP
        return dTdP

    def _TLherzLiquidus(self, P, **kwargs):
        """
        Returns the temperature of the lherzolite liquidus at any given pressure.
        Eqn(5). This is the temperature at which the rock would be completely
        molten if cpx was remained present for the entirety of melting.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        Tlzliq:   float
            Lherzolite liquidus temperature (degC).
        """
        TLherzLiquidus = (self.parameters['B1'] + self.parameters['B2']*P
                          + self.parameters['B3']*(P**2))
        return TLherzLiquidus

    def _RescaledTcpx(self, T, P, **kwargs):
        """
        Calculates the rescaled temperature during cpx-present melting. (Eq 3).

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
        TSolidus = self.TSolidus(P, **kwargs)
        TLherzLiquidus = self._TLherzLiquidus(P, **kwargs)
        RescaledTemperaturecpx = ((T - TSolidus)/(TLherzLiquidus - TSolidus))
        return RescaledTemperaturecpx

    def _Fcpx(self, T, P, **kwargs):
        """
        Melt fraction during cpx-present melting at the given P and T. Eq(2).

        Parameters
        ----------
        T : float
            Temperature (degC).
        P : float
            Pressure (degC).

        Returns
        -------
        float
            Melt fraction during cpx-present melting.
        """
        RescaledT = self._RescaledTcpx(T, P, **kwargs)
        Fcpx = RescaledT**self.parameters['beta1']
        return Fcpx

    def _RxnCoef(self, P, **kwargs):
        """
        Reaction coefficient for cpx during melting at the specified pressure. Eq(7).

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Reaction Coefficient.
        """
        RxnCoef = self.parameters['r1'] + self.parameters['r2']*P
        return RxnCoef

    def _FcpxOut(self, P, **kwargs):
        """
        Calculates the melt fraction required to exhaust cpx from the residue at the given
        pressure. Eq(6).

        Parameters
        ----------
        P : float
            Pressure (GPa)

        Returns
        -------
        float
            Melt fraction at which cpx is exhausted.
        """
        FcpxOut = self.parameters['Mcpx'] / self._RxnCoef(P, **kwargs)
        return FcpxOut

    def _TcpxOut(self, P, **kwargs):
        """
        Calculates the temperature at which cpx will be exhausted during melting at the given
        pressure. Eq(9).

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Temperature of cpx-exhaustion.
        """
        TSolidus = self.TSolidus(P, **kwargs)
        TcpxOut = (((self._FcpxOut(P, **kwargs)**(1/self.parameters['beta1'])))
                   * (self._TLherzLiquidus(P, **kwargs) - TSolidus) + TSolidus)
        return TcpxOut

    def _RescaledTopx(self, T, P, **kwargs):
        """
        Calculates the rescaled temperature during cpx-absent melting.

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
        TcpxOut = self._TcpxOut(P, **kwargs)
        RescaledTopx = ((T-TcpxOut) / (self.TLiquidus(P, **kwargs)-TcpxOut))
        return RescaledTopx

    def _Fopx(self, T, P, **kwargs):
        """
        Melt fraction during cpx-absent melting at the given P and T. Eq(8).

        Parameters
        ----------
        T : float
            Temperature (degC).
        P : float
            Pressure (degC).

        Returns
        -------
        float
            Melt fraction during cpx-absent melting.
        """
        FcpxOut = self._FcpxOut(P, **kwargs)
        FopxDry = (FcpxOut + (1-FcpxOut)
                   * self._RescaledTopx(T, P, **kwargs)**self.parameters['beta2'])
        return FopxDry

    def _dTdPSolidus(self, P, **kwargs):
        return self.parameters['A2'] + 2*self.parameters['A3']*P

    def _dTdPLiquidus(self, P, **kwargs):
        return self.parameters['C2'] + 2*self.parameters['C3']*P

    def _dTdPLherzLiquidus(self, P, **kwargs):
        return self.parameters['B2'] + 2*self.parameters['B3']*P

    def _dTdPcpxOut(self, P, **kwargs):
        term1 = ((self._TLherzLiquidus(P, **kwargs) - self.TSolidus(P, **kwargs))
                 / self.parameters['beta1']
                 * self._FcpxOut(P, **kwargs)**(1/self.parameters['beta1'] - 1))
        term2 = (self._FcpxOut(P, **kwargs)**(1/self.parameters['beta1'])
                 * (self._dTdPLherzLiquidus(P, **kwargs) - self._dTdPSolidus(P, **kwargs))
                 + self._dTdPSolidus(P, **kwargs))
        return term1 * self._dFdPcpxOut(P, **kwargs) + term2

    def _dFdPcpxOut(self, P, **kwargs):
        return - self.parameters['Mcpx'] / self._RxnCoef(P, **kwargs)**2 * self.parameters['r2']
