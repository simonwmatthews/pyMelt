"""
========================
Lambart et al. (Melt-PX)
========================

Implementation of the bulk-composition dependent Melt-PX
model by Lambart et al. (2016).
"""
from pyMelt.lithology_classes import lithology
from pyMelt.lithology_classes import default_properties
import numpy as np

class MeltPX(lithology):
    """
    The Melt-PX model by Lambart et al. (2016) for modelling
    a variety of pyroxenite bulk compositions (specified with
    parameters during initialisation).

    Parameters
    ----------
    alk : float
        The mole percent Na2O + K2O in the bulk composition
    Mgn : float
        The Mg# (Mg/(Mg + Fe)) on a mol basis, with all Fe as FeO
    Al2O3 : float
        The mole percent Al2O3 in the bulk composition
    CaO : float
        The mole percent CaO in the bulk composition
    high_F_productivity : float|None, default: None
        Extrapolate melt fraction expression beyond cpx-out (None)
        or used a fixed productivity (%/K).
    """
    def __init__(self,
                 alk : float,
                 Mgn : float,
                 Al2O3 : float,
                 CaO : float,
                 high_F_prod : float|None = None,
                 CP=default_properties['CP'],
                 alphas=default_properties['alphas'],
                 alphaf=default_properties['alphaf'],
                 rhos=default_properties['rhos'],
                 rhof=default_properties['rhof'],
                 DeltaS=default_properties['DeltaS'],
                 parameters={'a': 1.4,
                             'b': -7.08,
                             'c': -12.5,
                             'd': 150.0,
                             'e': 25.5,
                             'f': 1105.0,
                             'ap': -14.0,
                             'bp': 177,
                             'cp': 219,
                             'dp': -17.2,
                             'ep': 10.4,
                             'fp': 1000.0,
                             'Ap': -2.63,
                             'Bp': -1.58,
                             'Cp': 26.9,
                             'Dp': -7.5,
                             'Ep': -7.5,
                             'App': 2.86,
                             'Bpp': 18.0
                             },
                 ):
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.DeltaS = DeltaS
        self.parameters = parameters
        self.alk = alk
        self.Mgn = Mgn
        self.Al2O3 = Al2O3
        self.CaO = CaO
        self.high_F_prod = high_F_prod

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
        A = self._A(P)
        B = self._B(P)

        Tprime = (-B + np.sqrt(B**2 - 20 * A)) / 2 / A
        T_cpxout = self._T_cpxout(P)
        T_5percent = self._T_5percent(P)
        T = Tprime * (T_cpxout - T_5percent) + T_5percent

        return T

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
        A = self._A(P)
        B = self._B(P)

        Tprime = (-B + np.sqrt(B**2 + 4 * 95 * A)) / 2 / A
        T_cpxout = self._T_cpxout(P)
        T_5percent = self._T_5percent(P)
        T = Tprime * (T_cpxout - T_5percent) + T_5percent

        if T > T_cpxout and self.high_F_prod is not None:
            F_cpx_out = self._F_cpxout(P)
            dF = 1.0 - F_cpx_out
            dT = dF / self.high_F_prod
            T = T_cpxout + dT
        return T


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
        Tsol = self.TSolidus(P)
        Tliq = self.TLiquidus(P)
        assert Tliq > Tsol, "Out of model bounds, unphysical behaviour detected."
        if T < Tsol:
            return 0.0
        if T > Tliq:
            return 1.0
        
        T_cpxout = self._T_cpxout(P)
        T_5percent = self._T_5percent(P)

        if T > T_cpxout and self.high_F_prod is not None:
            F_cpx_out = self._F_cpxout(P)
            F = F_cpx_out + (T - T_cpxout) * self.high_F_prod
            return F
        

        Tprime = (T - T_5percent) / (T_cpxout - T_5percent)

        A = self._A(P)
        B = self._B(P)


        F = (A * Tprime**2 + B*Tprime + 5) / 100.0
        if F < 0.0:
            F = 0.0
        elif F > 1.0:
            F = 1.0
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
        Tsol = self.TSolidus(P)
        Tliq = self.TLiquidus(P)
        assert Tliq > Tsol, "Out of model bounds, unphysical behaviour detected."
        if T < Tsol or T > Tliq:
            return 0.0

        T_cpxout = self._T_cpxout(P)
        T_5percent = self._T_5percent(P)

        if T > T_cpxout and self.high_F_prod is not None:
            return 1.0 / self.high_F_prod

        Tprime = (T - T_5percent) / (T_cpxout - T_5percent)
        A = self._A(P) / 100.0
        B = self._B(P) / 100.0

        dTdF = (
            (T_cpxout - T_5percent) / (2 * A * Tprime + B)
        )

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
        Tsol = self.TSolidus(P)
        Tliq = self.TLiquidus(P)
        assert Tliq > Tsol, f"Out of model bounds, unphysical behaviour detected at P={P}."
        if T < Tsol:
            return self.alphas / self.rhos / self.CP
        elif T > Tliq:
            return self.alphaf / self.rhof / self.CP

        T_cpxout = self._T_cpxout(P)
        T_5percent = self._T_5percent(P)

        p = self.parameters
        Tprime = (T - T_5percent) / (T_cpxout - T_5percent)

        A = self._A(P) / 100.0
        B = self._B(P) / 100.0
        
        dAdP = (
            (2*(p['Ap'] + p['Bp'] * self.Mgn) 
             + p['Cp'] + p['Dp']*self.Mgn) / 100.0
        )

        dTcpxdP = 2 * p['ap'] * P + p['bp']

        if T > T_cpxout and self.high_F_prod is not None:
            dFcpxdP = (
                dAdP + 2 * A * Tprime
            )
            dTdP = dTcpxdP - dFcpxdP/self.high_F_prod
            
           # dTdP = np.nan
            return dTdP

       
        dTdP_prime = - Tprime**2 * dAdP / (2 * A * Tprime + B)
        dT5pdP = (
            2 * P * (p['a'] * self.alk + p['b'])
            + p['c'] * self.alk + p['d'] * self.Mgn + p['e']
        )
        dTdP = (
            (dTdP_prime 
             + (T - T_5percent)/(T_cpxout - T_5percent)**2 
             * (dTcpxdP - dT5pdP)
             + dT5pdP / (T_cpxout - T_5percent))
            * (T_cpxout - T_5percent)
        )
        return dTdP


        

    def _T_5percent(self, P):
        """
        Eqn 2.
        """
        p = self.parameters
        res = (
            P**2 * (p['a']*self.alk + p['b'])
            + P * (p['c']*self.alk + p['d']*self.Mgn + p['e'])
            + p['f']
        )

        return res

    def _T_cpxout(self, P):
        """
        Eqn. 3
        """
        p = self.parameters
        res = (
            p['ap']*P**2 + p['bp']*P + p['cp']*self.Mgn
            + p['dp']*self.Al2O3 + p['ep']*self.CaO + p['fp']
        )
        return res

    def _A(self, P):
        """
        Eqn. 5
        """
        p = self.parameters
        A = (
            (p['Ap'] + p['Bp'] * self.Mgn) * P**2
            + (p['Cp'] + p['Dp'] * self.Mgn) * P
            + p['Ep']
            )
        return A

    def _B(self, P):
        """
        Eqn. 6
        """
        p = self.parameters
        B = p['App']*self.CaO + p['Bpp']
        return B
    
    def _F_cpxout(self, P) -> float:
        """
        Calculate the melt fraction at which cpx is exhausted.

        Returns
        -------
        float
            Melt fraction.
        """
        Tsol = self.TSolidus(P)
        
        T_cpxout = self._T_cpxout(P)
        T_5percent = self._T_5percent(P)

        T = T_cpxout

        Tprime = (T - T_5percent) / (T_cpxout - T_5percent)

        A = self._A(P)
        B = self._B(P)


        F = (A * Tprime**2 + B*Tprime + 5) / 100.0

        return F

