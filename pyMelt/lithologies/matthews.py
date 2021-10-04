"""
======================
Matthews et al. (2021)
======================

Implentation of the melting models developed by Matthews et al. (2021).
"""

from pyMelt.lithology_classes import lithology as _lithology
from pyMelt.lithology_classes import default_properties as _default_properties


import numpy as _np


class kg1(_lithology):
    """
    Implementation of the KG1 melting model from Matthews et al. (2021).

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class.

    - Mcpx:   Mass fraction of cpx in the source. Controls the transition to
      low-productivity harzburgite-type melting.
    - A1:     Parameter used to define solidus.
    - A2:     Parameter used to define solidus.
    - A3:     Parameter used to define solidus.
    - A4:     Parameter used to define solidus.
    - B1:     Parameter used to define liquidus.
    - B2:     Parameter used to define liquidus.
    - B3:     Parameter used to define liquidus.
    - B4:     Parameter used to define liquidus.
    - C:      Parameter used to define lherzolite-liquidus.
    - beta1:  Parameter used to calculate melt fraction during cpx-present melting.
    - beta2:  Parameter used to calculate melt fraction during cpx-absent melting.
    - r1:     Parameter used to define cpx reaction coefficient.
    - r2:     Parameter used to define cpx reaction coefficient.

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    also be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: pyMelt.default_properties['CP']
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: pyMelt.default_properties['alphas']
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt.default_properties['alphaf']
        The thermal expansivity of the melt (1e-6 K-1)
    rhos :       float, default: pyMelt.default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt.default_properties['rhof']
        The density of the melt (kg m-3)
    DeltaS :     float, default: pyMelt.default_properties['DeltaS']
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
                 parameters={'A1':    450.000,
                             'A2':      2.098,
                             'A3':     17.000,
                             'A4':    623.828,
                             'B1':    174.566,
                             'B2':    336.833,
                             'B3':     66.762,
                             'B4':    503.101,
                             'C':       0.506,
                             'beta1':   1.382,
                             'beta2':   1.800,
                             'r1':      0.342,
                             'r2':      0.191,
                             'Mcpx':    0.342
                             }
                 ):
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.DeltaS = DeltaS
        self.parameters = parameters

    def TSolidus(self, P, **kwargs):
        """
        Returns the temperature of the solidus at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Solidus temperature (degC).
        """
        TSolidus = (self.parameters['A1']*_np.log(P + self.parameters['A2'])
                    + self.parameters['A3']*P
                    + self.parameters['A4'])
        return TSolidus

    def TLiquidus(self, P, **kwargs):
        """
        Returns the temperature of the liquidus at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Liquidus temperature (degC).
        """
        TLiquidus = (self.parameters['B1']*_np.log(P + self.parameters['B2'])
                     + self.parameters['B3']*P + self.parameters['B4'])
        return TLiquidus

    def F(self, P, T, **kwargs):
        """
        Wrapper for the melt fraction functions. If T and P are below the solidus,
        returns 0, if they are above the liquidus, returns 1. If below the temperature
        of cpx-exhaustion, calls the Fcpx function, otherwise calls the Fopx function.

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
        Calculates dT/dF(const. P). First calculates the melt fraction. If F is
        zero, returns _np.inf. If F is 1, returns _np.inf. Otherwise uses the
        appropriate expressions for cpx present or absent melting.

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
        P:  float
            Pressure (GPa).
        T:  float
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
        TcpxOut = self._dTdPcpxOut(P, **kwargs)
        dTdPcpxOut = self._dTdPcpxOut(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        dFdPcpxOut = self._dFdPcpxOut(P, **kwargs)

        if F == 0:
            dTdP = self.alphas/self.rhos/self.CP
        elif F < self._FcpxOut(P, **kwargs):
            dTdP = (((F**(1/self.parameters['beta1'])) * (dTdPLherzLiquidus-dTdPSolidus))
                    + dTdPSolidus)
        elif F < 1.0:
            Trel = (T - TcpxOut) / (TLiquidus - TcpxOut)
            dTdP = ((TLiquidus - TcpxOut)/(1-FcpxOut)
                    * (1/self.parameters['beta2'])*Trel**(1-self.parameters['beta2'])
                    * dFdPcpxOut * (Trel**self.parameters['beta2']-1)
                    + dTdPcpxOut + Trel*(dTdPLiquidus - dTdPcpxOut))
        else:
            dTdP = self.alphaf/self.rhof/self.CP
        return dTdP

    def _dTdPSolidus(self, P, **kwargs):
        """
        Returns the solidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Solidus temperaure gradient (degC/GPa)
        """
        dTdPSolidus = (self.parameters['A1']/(P + self.parameters['A2'])) + self.parameters['A3']
        return dTdPSolidus

    def _dTdPLiquidus(self, P, **kwargs):
        """
        Returns the liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Liquidus temperaure gradient (degC/GPa)
            """
        dTdPLiquidus = ((self.parameters['B1']/(P + self.parameters['B2']))
                        + self.parameters['B3'])
        return dTdPLiquidus

    def _TLherzLiquidus(self, P, **kwargs):
        """
        Returns the temperature of the lherzolite liquidus at any given pressure.
        This is the temperature at which the rock would be completely
        molten if cpx was remained present for the entirety of melting.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Lherzolite liquidus temperature (degC).
        """
        TSolidus = self.TSolidus(P, **kwargs)
        TLiquidus = self.TLiquidus(P, **kwargs)
        TLherzLiquidus = self.parameters['C']*TSolidus + (1 - self.parameters['C'])*TLiquidus
        return TLherzLiquidus

    def _dTdPLherzLiquidus(self, P, **kwargs):
        """
        Returns the lherzolite liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Lherzolite liquidus temperaure gradient (degC/GPa)
        """
        dTdPSolidus = self._dTdPSolidus(P, **kwargs)
        dTdPLiquidus = self._dTdPLiquidus(P, **kwargs)
        dTdPLherzoliteLiquidus = (self.parameters['C']*dTdPSolidus
                                  + (1 - self.parameters['C'])*dTdPLiquidus)
        return dTdPLherzoliteLiquidus

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
        RescaledTemperaturecpx = ((T - TSolidus)/(TLherzLiquidus-TSolidus))
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
        Calculates the melt fraction required to exhaust cpx from the residue
        at the given pressure. Eq(6).

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

    def _dFdPcpxOut(self, P, **kwargs):
        """
        Calculates the first derivative of FcpxOut.

        Parameters
        ----------
        P : float
            Pressure (GPa)

        Returns
        -------
        float
            First derivative of FcpxOut.
        """
        RxnCoef = self._RxnCoef(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        dFdPcpxOut = ((-1)*FcpxOut*self.parameters['r2'])/RxnCoef
        return dFdPcpxOut

    def _TcpxOut(self, P, **kwargs):
        """
        Calculates the temperature at which cpx will be exhausted during melting
        at the given pressure. Eq(9).

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
        TLherzLiquidus = self._TLherzLiquidus(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        TcpxOut = (((FcpxOut**(1/self.parameters['beta1']))) * (TLherzLiquidus - TSolidus)
                   + TSolidus)
        return TcpxOut

    def _dTdPcpxOut(self, P, **kwargs):
        """
        Calculates the temperature gradient of the cpx-exhaustion surface.

        Parameters
        ----------

        P : float
            Pressure (GPa)

        Returns
        -------
        float
            Temperature gradient of cpx-exhaustion.
        """
        TSolidus = self.TSolidus(P, **kwargs)
        TLherzLiquidus = self._TLherzLiquidus(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        dTdPSolidus = self._dTdPSolidus(P, **kwargs)
        dTdPLherzLiquidus = self._dTdPLherzLiquidus(P, **kwargs)
        dFdPcpxOut = self._dFdPcpxOut(P, **kwargs)
        A = ((dFdPcpxOut*(1/self.parameters['beta1'])
             * ((FcpxOut**((1/self.parameters['beta1'])-1)))) *
             (TLherzLiquidus-TSolidus))
        B = ((FcpxOut**(1/self.parameters['beta1']))*(dTdPLherzLiquidus-dTdPSolidus))+dTdPSolidus
        dTdPcpxOut = A + B
        return dTdPcpxOut

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
        TLiquidus = self.TLiquidus(P, **kwargs)
        RescaledTopx = ((T-TcpxOut)/(TLiquidus-TcpxOut))
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
        RescaledTopx = self._RescaledTopx(T, P, **kwargs)
        FopxDry = FcpxOut + (1-FcpxOut) * RescaledTopx**self.parameters['beta2']
        return FopxDry


class klb1:
    """
    Implementation of the KLB1 melting model from Matthews et al. (2021).

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class.

    - Mcpx:   Mass fraction of cpx in the source. Controls the transition to low-productivity
      harzburgite-type melting.
    - A1:     Parameter used to define solidus.
    - A2:     Parameter used to define solidus.
    - A3:     Parameter used to define solidus.
    - A4:		Parameter used to define solidus.
    - B1:     Parameter used to define liquidus.
    - B2:     Parameter used to define liquidus.
    - B3:     Parameter used to define liquidus.
    - B4:     Parameter used to define liquidus.
    - C:     	Parameter used to define lherzolite-liquidus.
    - beta1:  Parameter used to calculate melt fraction during cpx-present melting.
    - beta2:  Parameter used to calculate melt fraction during cpx-absent melting.
    - r1:     Parameter used to define cpx reaction coefficient.
    - r2:     Parameter used to define cpx reaction coefficient.

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    also be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: pyMelt.default_properties['CP']
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: pyMelt.default_properties['alphas']
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt.default_properties['alphaf']
        The thermal expansivity of the melt (1e-6 K-1)
    rhos :       float, default: pyMelt.default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt.default_properties['rhof']
        The density of the melt (kg m-3)
    DeltaS :     float, default: pyMelt.default_properties['DeltaS']
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
                 parameters={'Mcpx':     0.1500,
                             'A1':    2445.7540,
                             'A2':       9.5110,
                             'A3': - 99.7820,
                             'A4': - 4378.5810,
                             'B1':     480.4030,
                             'B2':     672.3910,
                             'B3':      12.2750,
                             'B4': - 1242.5360,
                             'C':        0.6873,
                             'beta1':    1.5000,
                             'beta2':    1.5000,
                             'r1':       0.5000,
                             'r2':       0.0800
                             }
                 ):
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.DeltaS = DeltaS
        self.parameters = parameters

    def TSolidus(self, P, **kwargs):
        """
        Returns the temperature of the solidus at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Solidus temperature (degC).
        """
        TSolidus = (self.parameters['A1']*_np.log(P + self.parameters['A2'])
                    + self.parameters['A3']*P + self.parameters['A4'])
        return TSolidus

    def TLiquidus(self, P, **kwargs):
        """
        Returns the temperature of the liquidus at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Liquidus temperature (degC).
        """
        TLiquidus = (self.parameters['B1']*_np.log(P + self.parameters['B2'])
                     + self.parameters['B3']*P + self.parameters['B4'])
        return TLiquidus

    def F(self, P, T, **kwargs):
        """
        Wrapper for the melt fraction functions. If T and P are below the solidus,
        returns 0, if they are above the liquidus, returns 1. If below the temperature
        of cpx-exhaustion, calls the Fcpx function, otherwise calls the Fopx function.

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
        Calculates dT/dF(const. P). First calculates the melt fraction. If F is
        zero, returns _np.inf. If F is 1, returns _np.inf. Otherwise uses the
        appropriate expressions for cpx present or absent melting.

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
                    * (self._TLherzLiquidus(P, **kwargs)-self.TSolidus(P, **kwargs))
                    * (F**((1-self.parameters['beta1'])/self.parameters['beta1'])))
        elif F < 1.0:
            dTdF = (((1/self.parameters['beta2']))
                    * (self.TLiquidus(P, **kwargs)-self._TcpxOut(P, **kwargs))
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
        TcpxOut = self._dTdPcpxOut(P, **kwargs)
        dTdPcpxOut = self._dTdPcpxOut(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        dFdPcpxOut = self._dFdPcpxOut(P, **kwargs)

        if F == 0:
            dTdP = self.alphas/self.rhos/self.CP
        elif F < self._FcpxOut(P, **kwargs):
            dTdP = (((F**(1/self.parameters['beta1']))
                    * (dTdPLherzLiquidus-dTdPSolidus)) + dTdPSolidus)
        elif F < 1.0:
            Trel = (T - TcpxOut) / (TLiquidus - TcpxOut)
            dTdP = ((TLiquidus - TcpxOut)/(1 - FcpxOut)
                    * (1/self.parameters['beta2'])*Trel**(1-self.parameters['beta2'])
                    * dFdPcpxOut * (Trel**self.parameters['beta2']-1)
                    + dTdPcpxOut + Trel*(dTdPLiquidus - dTdPcpxOut))
        else:
            dTdP = self.alphaf/self.rhof/self.CP
        return dTdP

    def _dTdPLiquidus(self, P, **kwargs):
        """
        Returns the liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Liquidus temperaure gradient (degC/GPa)
        """
        dTdPLiquidus = (self.parameters['B1']/(P + self.parameters['B2'])) + self.parameters['B3']
        return dTdPLiquidus

    def _TLherzLiquidus(self, P, **kwargs):
        """
        Returns the temperature of the lherzolite liquidus at any given pressure.
        This is the temperature at which the rock would be completely
        molten if cpx was remained present for the entirety of melting.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Lherzolite liquidus temperature (degC).
        """
        TSolidus = self.TSolidus(P, **kwargs)
        TLiquidus = self.TLiquidus(P, **kwargs)
        TLherzLiquidus = self.parameters['C']*TSolidus + (1 - self.parameters['C'])*TLiquidus
        return TLherzLiquidus

    def _dTdPSolidus(self, P, **kwargs):
        """
        Returns the solidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Solidus temperaure gradient (degC/GPa)
        """
        dTdPSolidus = (self.parameters['A1']/(P + self.parameters['A2'])) + self.parameters['A3']
        return dTdPSolidus

    def _dTdPLherzLiquidus(self, P, **kwargs):
        """
        Returns the lherzolite liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Lherzolite liquidus temperaure gradient (degC/GPa)
        """
        dTdPSolidus = self._dTdPSolidus(P, **kwargs)
        dTdPLiquidus = self._dTdPLiquidus(P, **kwargs)
        dTdPLherzoliteLiquidus = (self.parameters['C']*dTdPSolidus
                                  + (1 - self.parameters['C'])*dTdPLiquidus)
        return dTdPLherzoliteLiquidus

    def _RescaledTcpx(self, T, P, **kwargs):
        """
        Calculates the rescaled temperature during cpx-present melting.

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (GPa).

        Returns
        -------
        float
            Rescaled Temperature (dimensionless).
        """
        TSolidus = self.TSolidus(P, **kwargs)
        TLherzLiquidus = self._TLherzLiquidus(P, **kwargs)
        RescaledTemperaturecpx = ((T - TSolidus)/(TLherzLiquidus-TSolidus))
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
        Fcpx = self._RescaledTcpx(T, P, **kwargs)**self.parameters['beta1']
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
        Calculates the melt fraction required to exhaust cpx from the residue
        at the given pressure.

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

    def _dFdPcpxOut(self, P, **kwargs):
        """
        Calculates the first derivative of FcpxOut.

        Parameters
        ----------
        P : float
            Pressure (GPa)

        Returns
        -------
        float
            First derivative of FcpxOut.
        """
        RxnCoef = self._RxnCoef(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        return ((-1)*FcpxOut*self.parameters['r2'])/RxnCoef

    def _TcpxOut(self, P, **kwargs):
        """
        Calculates the temperature at which cpx will be exhausted during melting
        at the given pressure. Eq(9).

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
        TLherzLiquidus = self._TLherzLiquidus(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        TcpxOut = (((FcpxOut**(1/self.parameters['beta1'])))
                   * (TLherzLiquidus-TSolidus) + TSolidus)
        return TcpxOut

    def _dTdPcpxOut(self, P, **kwargs):
        """
        Calculates the temperature gradient of the cpx-exhaustion surface.

        Parameters
        ----------

        P : float
            Pressure (GPa)

        Returns
        -------
        float
            Temperature gradient of cpx-exhaustion.
        """
        TSolidus = self.TSolidus(P, **kwargs)
        TLherzLiquidus = self._TLherzLiquidus(P, **kwargs)
        FcpxOut = self._FcpxOut(P, **kwargs)
        dTdPSolidus = self._dTdPSolidus(P, **kwargs)
        dTdPLherzLiquidus = self._dTdPLherzLiquidus(P, **kwargs)
        dFdPcpxOut = self._dFdPcpxOut(P, **kwargs)

        A = ((dFdPcpxOut*(1/self.parameters['beta1'])
              * ((FcpxOut**((1/self.parameters['beta1'])-1))))
             * (TLherzLiquidus-TSolidus))
        B = (((FcpxOut**(1/self.parameters['beta1']))
             * (dTdPLherzLiquidus-dTdPSolidus))+dTdPSolidus)
        dTdPcpxOut = A + B
        return dTdPcpxOut

    def _RescaledTopx(self, T, P, **kwargs):
        """
        Calculates the rescaled temperature during cpx-absent melting.

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (GPa).

        Returns
        -------
        float
            Rescaled Temperature (dimensionless).
        """
        TcpxOut = self._TcpxOut(P, **kwargs)
        TLiquidus = self.TLiquidus(P, **kwargs)
        RescaledTopx = ((T-TcpxOut)/(TLiquidus-TcpxOut))
        return RescaledTopx

    def _Fopx(self, T, P, **kwargs):
        """
        Melt fraction during cpx-absent melting at the given P and T.

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (degC).

        Returns
        -------
        float
            Melt fraction during cpx-absent melting.
        """
        FcpxOut = self._FcpxOut(P, **kwargs)
        RescaledTopx = self._RescaledTopx(T, P, **kwargs)
        FopxDry = FcpxOut + (1-FcpxOut) * RescaledTopx**self.parameters['beta2']
        return FopxDry


class eclogite(_lithology):
    """
    Implementation of the silica-saturated pyroxenite (or eclogite) melting model from
    Matthews et al. (2021).

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class.

    - C1:     Parameter used in solidus definition.
    - C2:     Parameter used in solidus definition.
    - C3:     Parameter used in solidus definition.
    - C4:     Parameter used in solidus definition.
    - D1:     Parameter used in liquidus definition.
    - D2:     Parameter used in liquidus definition.
    - D3:     Parameter used in liquidus definition.
    - D4:     Parameter used in liquidus definition.
    - beta:   Parameter used in melt fraction definition.

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    also be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: pyMelt.default_properties['CP']
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: pyMelt.default_properties['alphas']
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt.default_properties['alphaf']
        The thermal expansivity of the melt (1e-6 K-1)
    rhos :       float, default: pyMelt.default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt.default_properties['rhof']
        The density of the melt (kg m-3)
    DeltaS :     float, default: pyMelt.default_properties['DeltaS']
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
                 parameters={'C1':    533.842,
                             'C2':      4.921,
                             'C3':     20.148,
                             'C4':     80.879,
                             'D1':    994.149,
                             'D2':      8.092,
                             'D3': - 11.778,
                             'D4': - 862.641,
                             'beta':    2.134,
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
        Calculates melt fraction at a given pressure and temperature using T'**beta, where T is the
        normalised temperature: (T-Tsolidus)/(T-Tliquidus). If P and T are below the the solidus,
        0 is returned, if they are above the liquidus, 1 is returned.

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
            F = Tr**self.parameters['beta']
        return F

    def TSolidus(self, P, **kwargs):
        """
        Calculates the solidus temperature at a given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Solidus temperature (degC).
        """
        Tsol = (self.parameters['C1']
                * _np.log(P + self.parameters['C2'])
                + self.parameters['C3']*P + self.parameters['C4'])
        return Tsol

    def TLiquidus(self, P, **kwargs):
        """
        Calculates the liquidus temperature at a given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa)

        Returns
        -------
        float
            Liquidus temperature (degC).
        """
        Tliq = (self.parameters['D1'] * _np.log(P + self.parameters['D2'])
                + self.parameters['D3']*P + self.parameters['D4'])
        return Tliq

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
        F = self.F(P, T, **kwargs)
        if T < Tsol:
            dTdF = _np.inf
        elif T > Tliq:
            dTdF = _np.inf
        else:
            dTdF = ((1/self.parameters['beta']) * (Tliq-Tsol) * F**((1/self.parameters['beta'])-1))

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
        dTdPsol = self._dTdPSolidus(P, **kwargs)
        dTdPliq = self._dTdPLiquidus(P, **kwargs)
        F = self.F(P, T, **kwargs)
        if F == 0:
            dTdP = self.alphas/self.rhos/self.CP
        elif F == 1:
            dTdP = self.alphaf/self.rhof/self.CP
        else:
            dTdP = (F**(1/self.parameters['beta'])) * (dTdPliq-dTdPsol) + dTdPsol
        return dTdP

    def _dTdPSolidus(self, P, **kwargs):
        """
        Returns the solidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Solidus temperaure gradient (degC/GPa)
        """
        dTdPSolidus = ((self.parameters['C1']/(P + self.parameters['C2']))
                       + self.parameters['C3'])
        return dTdPSolidus

    def _dTdPLiquidus(self, P, **kwargs):
        """
        Returns the liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P : float
            Pressure (GPa).

        Returns
        -------
        float
            Liquidus temperaure gradient (degC/GPa)
        """
        dTdPLiquidus = ((self.parameters['D1']/(P + self.parameters['D2']))
                        + self.parameters['D3'])
        return dTdPLiquidus
