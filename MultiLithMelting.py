#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 20:25:26 2017

@author: simon
"""
# Version 1.4. 4 November 2019. Updated the new melting parameterisations to exclude
#                               transition zone melting and to incorporate more
#                               pyroxenite melting experiments.
# Version 1.3. 2 March 2019. Incorporated isobaric melting for super-solidus starts.
#                            Incorporated the new parameterisations for KLB1, G2 and KG1.
# Version 1.2. 13 November 2018. Fixed bugs occuring some of the time when lithologies
#                                are exhausted.
# Version 1.1. 07 September 2018. Added extracted melt temperature calculation.

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
#params = {'text.usetex': False, 'mathtext.fontset': 'stixsans',
#          'xtick.labelsize':10, 'ytick.labelsize':10, 'axes.labelsize':10,
#          'xtick.direction':'in','ytick.direction':'in','figure.figsize':(6,4),
#          'xtick.top':True,'ytick.right':True}
#plt.rcParams.update(params)



class LithologyKG1:
    """
    Lithology formatted like our new lherzolite melting model. Default parameter
    values taken from Katz et al. (2003) for dry lherzolite melting. All parameters
    listed below can be called. Calibrated on thermocalc results for KG1 up to
    70 kbar.

    Parameters
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
    Mcpx:   float
        Mass fraction of cpx in the source. Controls the transition to
        low-productivity harzburgite-type melting.
    A1:     float
        Parameter used to define solidus.
    A2:     float
        Parameter used to define solidus.
    A3:     float
        Parameter used to define solidus.
	A4:		float
		Parameter used to define solidus.
    B1:     float
        Parameter used to define liquidus.
    B2:     float
        Parameter used to define liquidus.
    B3:     float
        Parameter used to define liquidus.
	B4:     float
        Parameter used to define liquidus.
    C:     	float
        Parameter used to define lherzolite-liquidus.
    beta1:  float
        Parameter used to calculate melt fraction during cpx-present melting.
    beta2:  float
        Parameter used to calculate melt fraction during cpx-absent melting.
    r1:     float
        Parameter used to define cpx reaction coefficient.
    r2:     float
        Parameter used to define cpx reaction coefficient.
    """
#    def __init__(self,Mcpx=0.342,DeltaS=300.0,A1=519.458,A2=2.098,A3=-11.365,A4=623.828,B1=174.566,
#                 B2=336.833,B3=66.762,B4=503.101,C=0.506,beta1=1.382,beta2=1.371,
#                 r1=0.342,r2=0.191,CP=1000.0,alphas=40.0,alphaf=68.0,rhos=3.3,rhof=2.9):
#    x = [0.342,0.506,0.342,0.191,1.382,1.8,450,2.098,17,623.8,174.566,336.833,66.762,503.101]

    def __init__(self,Mcpx=0.342,DeltaS=300.0,A1=450,A2=2.098,A3=17,A4=623.828,B1=174.566,
                 B2=336.833,B3=66.762,B4=503.101,C=0.506,beta1=1.382,beta2=1.8,
                 r1=0.342,r2=0.191,CP=1000.0,alphas=40.0,alphaf=68.0,rhos=3.3,rhof=2.9):

        self.Mcpx = Mcpx
        self.DeltaS = DeltaS
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.A4 = A4
        self.B1 = B1
        self.B2 = B2
        self.B3 = B3
        self.B4 = B4
        self.C = C
        self.beta1 = beta1
        self.beta2 = beta2
        self.r1 = r1
        self.r2 = r2
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof

    def TSolidus(self,P):
        """
        Returns the temperature of the solidus at any given pressure.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tsol:   float
            Solidus temperature (degC).
        """
        _TSolidus = self.A1*np.log(P + self.A2) + self.A3*P + self.A4
        return _TSolidus

    def dTdPSolidus(self,P):
        """
        Returns the solidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPsol:	float
            Solidus temperaure gradient (degC/GPa)
        """
        _dTdPSolidus = (self.A1/(P + self.A2)) + self.A3
        return _dTdPSolidus

    def TLiquidus(self,P):
        """
        Returns the temperature of the liquidus at any given pressure.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tliq:   float
            Liquidus temperature (degC).
        """
        _TLiquidus = self.B1*np.log(P + self.B2) + self.B3*P + self.B4
        return _TLiquidus

    def dTdPLiquidus(self,P):
        """
        Returns the liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPliq:	float
            Liquidus temperaure gradient (degC/GPa)
            """
        _dTdPLiquidus = (self.B1/(P + self.B2)) + self.B3
        return _dTdPLiquidus

    def TLherzLiquidus(self,P):
        """
        Returns the temperature of the lherzolite liquidus at any given pressure.
        This is the temperature at which the rock would be completely
        molten if cpx was remained present for the entirety of melting.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tlzliq:   float
            Lherzolite liquidus temperature (degC).
        """
        _TSolidus = self.TSolidus(P)
        _TLiquidus = self.TLiquidus(P)
        _TLherzLiquidus = self.C*_TSolidus + (1 - self.C)*_TLiquidus
        return _TLherzLiquidus

    def dTdPLherzLiquidus(self,P):
        """
        Returns the lherzolite liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPlzliq:	float
            Lherzolite liquidus temperaure gradient (degC/GPa)
        """
        _dTdPSolidus = self.dTdPSolidus(P)
        _dTdPLiquidus = self.dTdPLiquidus(P)
        _dTdPLherzoliteLiquidus = self.C*_dTdPSolidus + (1 - self.C)*_dTdPLiquidus
        return _dTdPLherzoliteLiquidus

    def RescaledTcpx(self,T,P):
        """
        Calculates the rescaled temperature during cpx-present melting. (Eq 3).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (GPa).

        Returns
        -------
        T':  float
            Rescaled Temperature (dimensionless).
        """
        _TSolidus = self.TSolidus(P)
        _TLherzLiquidus = self.TLherzLiquidus(P)
        _RescaledTemperaturecpx = ((T - _TSolidus)/(_TLherzLiquidus-_TSolidus))
        return _RescaledTemperaturecpx

    def Fcpx(self,T,P):
        """
        Melt fraction during cpx-present melting at the given P and T. Eq(2).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (degC).

        Returns
        -------
        F:  float
            Melt fraction during cpx-present melting.
        """
        _RescaledT = self.RescaledTcpx(T,P)
        _Fcpx = _RescaledT**self.beta1
        return _Fcpx

    def RxnCoef(self,P):
        """
        Reaction coefficient for cpx during melting at the specified pressure. Eq(7).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        RxnCoef:    float
            Reaction Coefficient.
        """
        _RxnCoef = self.r1 + self.r2*P
        return _RxnCoef

    def FcpxOut(self,P):
        """
        Calculates the melt fraction required to exhaust cpx from the residue
        at the given pressure. Eq(6).

        Parameters
        ----------
        P:  float
            Pressure (GPa)

        Returns
        -------
        Fcpx-out:   float
            Melt fraction at which cpx is exhausted.
        """
        _FcpxOut = self.Mcpx / self.RxnCoef(P)
        return _FcpxOut

    def dFdPcpxOut(self,P):
        """
        Calculates the first derivative of FcpxOut.

        Parameters
        ----------
        P:	float
            Pressure (GPa)

        Returns
        -------
        dFdPcpx-out:	float
            First derivative of FcpxOut.
        """
        _RxnCoef = self.RxnCoef(P)
        _FcpxOut = self.FcpxOut(P)
        _dFdPcpxOut = ((-1)*_FcpxOut*self.r2)/_RxnCoef
        return _dFdPcpxOut

    def TcpxOut(self,P):
        """
        Calculates the temperature at which cpx will be exhausted during melting
        at the given pressure. Eq(9).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tcpx-out:   float
            Temperature of cpx-exhaustion.
        """
        _TSolidus = self.TSolidus(P)
        _TLherzLiquidus = self.TLherzLiquidus(P)
        _FcpxOut = self.FcpxOut(P)
        _TcpxOut= ((_FcpxOut**(1/self.beta1)))*(_TLherzLiquidus-_TSolidus)+_TSolidus
        return _TcpxOut

    def dTdPcpxOut(self,P):
        """
        Calculates the temperature gradient of the cpx-exhaustion surface.

        Parameters
        ----------

        P:	float
            Pressure (GPa)

        Returns
        -------
        dTdPcpx-out:	float
            Temperature gradient of cpx-exhaustion.
        """
        _TSolidus = self.TSolidus(P)
        _TLherzLiquidus = self.TLherzLiquidus(P)
        _FcpxOut = self.FcpxOut(P)
        _dTdPSolidus = self.dTdPSolidus(P)
        _dTdPLherzLiquidus = self.dTdPLherzLiquidus(P)
        _dFdPcpxOut = self.dFdPcpxOut(P)
        _A = (_dFdPcpxOut*(1/self.beta1)*((_FcpxOut**((1/self.beta1)-1))))*(_TLherzLiquidus-_TSolidus)
        _B = ((_FcpxOut**(1/self.beta1))*(_dTdPLherzLiquidus-_dTdPSolidus))+_dTdPSolidus
        _dTdPcpxOut = _A + _B
        return _dTdPcpxOut

    def RescaledTopx(self,T,P):
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
        T':  float
            Rescaled Temperature (dimensionless).
        """
        _TcpxOut = self.TcpxOut(P)
        _TLiquidus = self.TLiquidus(P)
        _RescaledTopx = ((T-_TcpxOut)/(_TLiquidus-_TcpxOut))
        return _RescaledTopx

    def Fopx(self,T,P):
        """
        Melt fraction during cpx-absent melting at the given P and T. Eq(8).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (degC).

        Returns
        -------
        F:  float
            Melt fraction during cpx-absent melting.
        """
        _FcpxOut = self.FcpxOut(P)
        _RescaledTopx = self.RescaledTopx(T,P)
        _FopxDry = _FcpxOut + (1-_FcpxOut)* _RescaledTopx**self.beta2
        return _FopxDry

    def F(self,P,T):
        """
        Wrapper for the melt fraction functions. If T and P are below the solidus,
        returns 0, if they are above the liquidus, returns 1. If below the temperature
        of cpx-exhaustion, calls the Fcpx function, otherwise calls the Fopx function.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        F:  float
            Melt fraction.
        """
        if T > self.TLiquidus(P):
            _F = 1.0
        elif T < self.TSolidus(P):
            _F = 0.0
        elif T < self.TcpxOut(P):
            _F = self.Fcpx(T,P)
        else:
            _F = self.Fopx(T,P)

        return _F

    def dTdF(self,P,T):
        """
        Calculates dT/dF(const. P). First calculates the melt fraction. If F is
        zero, returns np.inf. If F is 1, returns np.inf. Otherwise uses the
        appropriate expressions for cpx present or absent melting.

        Parameters
        ----------
        P:  float
            Pressure (GPa)
        T:  float
            Temperature (degC)

        Returns
        -------
        dTdF:   float
            dT/dF(const. P) (K).
        """
        _F = self.F(P,T)
        if _F == 0:
            dTdF = np.inf # If no melt fraction the derivative is zero. Prevents division by zero.
        elif _F < self.FcpxOut(P):
            dTdF = ((1/self.beta1))*(self.TLherzLiquidus(P)-self.TSolidus(P))*(_F**((1-self.beta1)/self.beta1))
        elif _F < 1.0:
            dTdF = ((1/self.beta2))*(self.TLiquidus(P)-self.TcpxOut(P))*(_F**((1-self.beta2)/self.beta2))
        else:
            dTdF = np.inf
        return dTdF


    def dTdP(self,P,T):
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
        dTdP:   float
            dT/dP(const. F) (K GPa-1).
        """
        _F = self.F(P,T)
        _FcpxOut = self.FcpxOut(P)
        _dTdPSolidus = self.dTdPSolidus(P)
        _TLiquidus = self.TLiquidus(P)
        _dTdPLiquidus = self.dTdPLiquidus(P)
        _dTdPLherzLiquidus = self.dTdPLherzLiquidus(P)
        _TcpxOut = self.dTdPcpxOut(P)
        _dTdPcpxOut = self.dTdPcpxOut(P)
        _FcpxOut = self.FcpxOut(P)
        _dFdPcpxOut = self.dFdPcpxOut(P)

        if _F == 0:
            _dTdP = self.alphas/self.rhos/self.CP
        elif _F < self.FcpxOut(P):
            _dTdP = ((_F**(1/self.beta1))*(_dTdPLherzLiquidus-_dTdPSolidus)) + _dTdPSolidus
        elif _F < 1.0:
            _Trel = (T- _TcpxOut)/(_TLiquidus-_TcpxOut)
            _dTdP = (_TLiquidus - _TcpxOut)/(1-_FcpxOut) * (1/self.beta2)*_Trel**(1-self.beta2) * _dFdPcpxOut * (_Trel**self.beta2-1) \
                    + _dTdPcpxOut + _Trel*(_dTdPLiquidus - _dTdPcpxOut)
        else:
            _dTdP = self.alphaf/self.rhof/self.CP
        return _dTdP





class LithologyKLB1:
    """
    Lithology formatted like our new lherzolite melting model. Does
    not incorporate their hydrous melting parameterisation. Default parameter
    values taken from Katz et al. (2003) for dry lherzolite melting. All parameters
    listed below can be called.

    Parameters
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
    Mcpx:   float
        Mass fraction of cpx in the source. Controls the transition to
        low-productivity harzburgite-type melting.
    A1:     float
        Parameter used to define solidus.
    A2:     float
        Parameter used to define solidus.
    A3:     float
        Parameter used to define solidus.
	A4:		float
		Parameter used to define solidus.
    B1:     float
        Parameter used to define liquidus.
    B2:     float
        Parameter used to define liquidus.
    B3:     float
        Parameter used to define liquidus.
	B4:     float
        Parameter used to define liquidus.
    C:     	float
        Parameter used to define lherzolite-liquidus.
    beta1:  float
        Parameter used to calculate melt fraction during cpx-present melting.
    beta2:  float
        Parameter used to calculate melt fraction during cpx-absent melting.
    r1:     float
        Parameter used to define cpx reaction coefficient.
    r2:     float
        Parameter used to define cpx reaction coefficient.
    """
    def __init__(self,Mcpx=0.15,DeltaS=300.0,A1=2445.754,A2=9.511,A3=-99.782,A4=-4378.581,B1=480.403,
                 B2=672.391,B3=12.275,B4=-1242.536,C=0.6873,beta1=1.5,beta2=1.5,
                 r1=0.5,r2=0.08,CP=1000.0,alphas=40.0,alphaf=68.0,rhos=3.3,rhof=2.9):
        self.Mcpx = Mcpx
        self.DeltaS = DeltaS
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.A4 = A4
        self.B1 = B1
        self.B2 = B2
        self.B3 = B3
        self.B4 = B4
        self.C = C
        self.beta1 = beta1
        self.beta2 = beta2
        self.r1 = r1
        self.r2 = r2
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof

    def TSolidus(self,P):
        """
        Returns the temperature of the solidus at any given pressure.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tsol:   float
            Solidus temperature (degC).
        """
        _TSolidus = self.A1*np.log(P + self.A2) + self.A3*P + self.A4
        return _TSolidus

    def dTdPSolidus(self,P):
        """
        Returns the solidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPsol:	float
            Solidus temperaure gradient (degC/GPa)
        """
        _dTdPSolidus = (self.A1/(P + self.A2)) + self.A3
        return _dTdPSolidus

    def TLiquidus(self,P):
        """
        Returns the temperature of the liquidus at any given pressure.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tliq:   float
            Liquidus temperature (degC).
        """
        _TLiquidus = self.B1*np.log(P + self.B2) + self.B3*P + self.B4
        return _TLiquidus

    def dTdPLiquidus(self,P):
        """
        Returns the liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPliq:	float
            Liquidus temperaure gradient (degC/GPa)
            """
        _dTdPLiquidus = (self.B1/(P + self.B2)) + self.B3
        return _dTdPLiquidus

    def TLherzLiquidus(self,P):
        """
        Returns the temperature of the lherzolite liquidus at any given pressure.
        This is the temperature at which the rock would be completely
        molten if cpx was remained present for the entirety of melting.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tlzliq:   float
            Lherzolite liquidus temperature (degC).
        """
        _TSolidus = self.TSolidus(P)
        _TLiquidus = self.TLiquidus(P)
        _TLherzLiquidus = self.C*_TSolidus + (1 - self.C)*_TLiquidus
        return _TLherzLiquidus

    def dTdPLherzLiquidus(self,P):
        """
        Returns the lherzolite liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPlzliq:	float
            Lherzolite liquidus temperaure gradient (degC/GPa)
        """
        _dTdPSolidus = self.dTdPSolidus(P)
        _dTdPLiquidus = self.dTdPLiquidus(P)
        _dTdPLherzoliteLiquidus = self.C*_dTdPSolidus + (1 - self.C)*_dTdPLiquidus
        return _dTdPLherzoliteLiquidus

    def RescaledTcpx(self,T,P):
        """
        Calculates the rescaled temperature during cpx-present melting. (Eq 3).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (GPa).

        Returns
        -------
        T':  float
            Rescaled Temperature (dimensionless).
        """
        _TSolidus = self.TSolidus(P)
        _TLherzLiquidus = self.TLherzLiquidus(P)
        _RescaledTemperaturecpx = ((T - _TSolidus)/(_TLherzLiquidus-_TSolidus))
        return _RescaledTemperaturecpx

    def Fcpx(self,T,P):
        """
        Melt fraction during cpx-present melting at the given P and T. Eq(2).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (degC).

        Returns
        -------
        F:  float
            Melt fraction during cpx-present melting.
        """
        _RescaledT = self.RescaledTcpx(T,P)
        _Fcpx = _RescaledT**self.beta1
        return _Fcpx

    def RxnCoef(self,P):
        """
        Reaction coefficient for cpx during melting at the specified pressure. Eq(7).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        RxnCoef:    float
            Reaction Coefficient.
        """
        _RxnCoef = self.r1 + self.r2*P
        return _RxnCoef

    def FcpxOut(self,P):
        """
        Calculates the melt fraction required to exhaust cpx from the residue
        at the given pressure. Eq(6).

        Parameters
        ----------
        P:  float
            Pressure (GPa)

        Returns
        -------
        Fcpx-out:   float
            Melt fraction at which cpx is exhausted.
        """
        _FcpxOut = self.Mcpx / self.RxnCoef(P)
        return _FcpxOut

    def dFdPcpxOut(self,P):
        """
        Calculates the first derivative of FcpxOut.

        Parameters
        ----------
        P:	float
            Pressure (GPa)

        Returns
        -------
        dFdPcpx-out:	float
            First derivative of FcpxOut.
        """
        _RxnCoef = self.RxnCoef(P)
        _FcpxOut = self.FcpxOut(P)
        _dFdPcpxOut = ((-1)*_FcpxOut*self.r2)/_RxnCoef
        return _dFdPcpxOut

    def TcpxOut(self,P):
        """
        Calculates the temperature at which cpx will be exhausted during melting
        at the given pressure. Eq(9).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tcpx-out:   float
            Temperature of cpx-exhaustion.
        """
        _TSolidus = self.TSolidus(P)
        _TLherzLiquidus = self.TLherzLiquidus(P)
        _FcpxOut = self.FcpxOut(P)
        _TcpxOut= ((_FcpxOut**(1/self.beta1)))*(_TLherzLiquidus-_TSolidus)+_TSolidus
        return _TcpxOut

    def dTdPcpxOut(self,P):
        """
        Calculates the temperature gradient of the cpx-exhaustion surface.

        Parameters
        ----------

        P:	float
            Pressure (GPa)

        Returns
        -------
        dTdPcpx-out:	float
            Temperature gradient of cpx-exhaustion.
        """
        _TSolidus = self.TSolidus(P)
        _TLherzLiquidus = self.TLherzLiquidus(P)
        _FcpxOut = self.FcpxOut(P)
        _dTdPSolidus = self.dTdPSolidus(P)
        _dTdPLherzLiquidus = self.dTdPLherzLiquidus(P)
        _dFdPcpxOut = self.dFdPcpxOut(P)
        _A = (_dFdPcpxOut*(1/self.beta1)*((_FcpxOut**((1/self.beta1)-1))))*(_TLherzLiquidus-_TSolidus)
        _B = ((_FcpxOut**(1/self.beta1))*(_dTdPLherzLiquidus-_dTdPSolidus))+_dTdPSolidus
        _dTdPcpxOut = _A + _B
        return _dTdPcpxOut

    def RescaledTopx(self,T,P):
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
        T':  float
            Rescaled Temperature (dimensionless).
        """
        _TcpxOut = self.TcpxOut(P)
        _TLiquidus = self.TLiquidus(P)
        _RescaledTopx = ((T-_TcpxOut)/(_TLiquidus-_TcpxOut))
        return _RescaledTopx

    def Fopx(self,T,P):
        """
        Melt fraction during cpx-absent melting at the given P and T. Eq(8).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (degC).

        Returns
        -------
        F:  float
            Melt fraction during cpx-absent melting.
        """
        _FcpxOut = self.FcpxOut(P)
        _RescaledTopx = self.RescaledTopx(T,P)
        _FopxDry = _FcpxOut + (1-_FcpxOut)* _RescaledTopx**self.beta2
        return _FopxDry

    def F(self,P,T):
        """
        Wrapper for the melt fraction functions. If T and P are below the solidus,
        returns 0, if they are above the liquidus, returns 1. If below the temperature
        of cpx-exhaustion, calls the Fcpx function, otherwise calls the Fopx function.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        F:  float
            Melt fraction.
        """
        if T > self.TLiquidus(P):
            _F = 1.0
        elif T < self.TSolidus(P):
            _F = 0.0
        elif T < self.TcpxOut(P):
            _F = self.Fcpx(T,P)
        else:
            _F = self.Fopx(T,P)

        return _F

    def dTdF(self,P,T):
        """
        Calculates dT/dF(const. P). First calculates the melt fraction. If F is
        zero, returns np.inf. If F is 1, returns np.inf. Otherwise uses the
        appropriate expressions for cpx present or absent melting.

        Parameters
        ----------
        P:  float
            Pressure (GPa)
        T:  float
            Temperature (degC)

        Returns
        -------
        dTdF:   float
            dT/dF(const. P) (K).
        """
        _F = self.F(P,T)
        if _F == 0:
            dTdF = np.inf # If no melt fraction the derivative is zero. Prevents division by zero.
        elif _F < self.FcpxOut(P):
            dTdF = ((1/self.beta1))*(self.TLherzLiquidus(P)-self.TSolidus(P))*(_F**((1-self.beta1)/self.beta1))
        elif _F < 1.0:
            dTdF = ((1/self.beta2))*(self.TLiquidus(P)-self.TcpxOut(P))*(_F**((1-self.beta2)/self.beta2))
        else:
            dTdF = np.inf
        return dTdF


    def dTdP(self,P,T):
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
        dTdP:   float
            dT/dP(const. F) (K GPa-1).
        """
        _F = self.F(P,T)
        _FcpxOut = self.FcpxOut(P)
        _dTdPSolidus = self.dTdPSolidus(P)
        _TLiquidus = self.TLiquidus(P)
        _dTdPLiquidus = self.dTdPLiquidus(P)
        _dTdPLherzLiquidus = self.dTdPLherzLiquidus(P)
        _TcpxOut = self.dTdPcpxOut(P)
        _dTdPcpxOut = self.dTdPcpxOut(P)
        _FcpxOut = self.FcpxOut(P)
        _dFdPcpxOut = self.dFdPcpxOut(P)

        if _F == 0:
            _dTdP = self.alphas/self.rhos/self.CP
        elif _F < self.FcpxOut(P):
            _dTdP = ((_F**(1/self.beta1))*(_dTdPLherzLiquidus-_dTdPSolidus)) + _dTdPSolidus
        elif _F < 1.0:
            _Trel = (T- _TcpxOut)/(_TLiquidus-_TcpxOut)
            _dTdP = (_TLiquidus - _TcpxOut)/(1-_FcpxOut) * (1/self.beta2)*_Trel**(1-self.beta2) * _dFdPcpxOut * (_Trel**self.beta2-1) \
                    + _dTdPcpxOut + _Trel*(_dTdPLiquidus - _dTdPcpxOut)
        else:
            _dTdP = self.alphaf/self.rhof/self.CP
        return _dTdP

class LithologyPx:
    """
    Lithology formatted like our silica-saturated pyroxenite melting model.
    Default values reproduce the behaviour of silica-saturated pyroxenite. Default thermodynamic
    constants are those used by Katz et al. (2003). All the parameters listed below
    are callable.

    Parameters
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
    C1:  float
        Parameter used in solidus definition.
    C2:  float
        Parameter used in solidus definition.
    C3:  float
        Parameter used in solidus definition.
    C4:  float
        Parameter used in solidus definition.
    D1:  float
        Parameter used in liquidus definition.
    D2:  float
        Parameter used in liquidus definition.
    D3:  float
        Parameter used in liquidus definition.
    D4:  float
        Parameter used in liquidus definition.
    beta:   float
        Parameter used in melt fraction definition.

    """
    def __init__(self,C1=533.842,C2=4.921,C3=20.148,C4=80.879,D1=994.149,
                 D2=8.092,D3=-11.778,D4=-862.641,DeltaS=300.0,CP=1000.0,alphas=40.0,alphaf=68.0,rhos=3.3,rhof=2.9,beta=2.134):
        self.C1=C1
        self.C2=C2
        self.C3=C3
        self.C4=C4
        self.D1=D1
        self.D2=D2
        self.D3=D3
        self.D4=D4
        self.DeltaS=DeltaS
        self.CP=CP
        self.alphas=alphas
        self.alphaf=alphaf
        self.rhos=rhos
        self.rhof=rhof
        self.beta=beta

    def F(self,P,T):
        """
        Calculates melt fraction at a given pressure and temperature using
        T'**beta, where T is the
        normalised temperature: (T-Tsolidus)/(T-Tliquidus). If P and T are
        below the the solidus, 0 is returned, if they are above the liquidus,
        1 is returned.

        Parameters
        ----------
        P:  float
            Pressure (GPa)
        T:  float
            Temperature (degC)

        Returns
        -------
        F:  float
            Melt fraction.
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        if T < _Tsol:
            _F = 0.0
        elif T > _Tliq:
            _F = 1.0
        else:
            _Tr = (T-_Tsol)/(_Tliq-_Tsol)
            _F = _Tr**self.beta
        return _F

    def TSolidus(self,P):
        """
        Calculates the solidus temperature at a given pressure.

        Parameters
        ----------
        P: float
            Pressure (GPa).

        Returns
        -------
        Tsol:   float
            Solidus temperature (degC).
        """
        _Tsol = self.C1*np.log(P + self.C2) + self.C3*P + self.C4
        return _Tsol

    def dTdPSolidus(self,P):
        """
        Returns the solidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPsol:	float
            Solidus temperaure gradient (degC/GPa)
        """
        _dTdPSolidus = (self.C1/(P + self.C2)) + self.C3
        return _dTdPSolidus

    def TLiquidus(self,P):
        """
        Calculates the liquidus temperature at a given pressure.

        Parameters
        ----------
        P:  float
            Pressure (GPa)

        Returns
        -------
        Tliq:   float
            Liquidus temperature (degC).
        """
        _Tliq = self.D1*np.log(P + self.D2) + self.D3*P + self.D4
        return(_Tliq)

    def dTdPLiquidus(self,P):
        """
        Returns the liquidus temperature gradient at any given pressure.

        Parameters
        ----------
        P:	float
            Pressure (GPa).

        Returns
        -------
        dTdPliq:	float
            Liquidus temperaure gradient (degC/GPa)
        """
        _dTdPLiquidus = (self.D1/(P + self.D2)) + self.D3
        return _dTdPLiquidus

    def dTdF(self,P,T):
        """
        Calculates dT/dF(const. P) at a given pressure and temperature.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        dTdF:   float
            dT/dF(const. P) (K)
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        _F = self.F(P,T)
        if T < _Tsol:
            _dTdF = np.inf
        elif T > _Tliq:
            _dTdF = np.inf
        else:
            _dTdF = (1/self.beta)*(_Tliq-_Tsol)*_F**((1/self.beta)-1)

        return _dTdF

    def dTdP(self,P,T):
        """
        Calculates dT/dP(const. F) at a given pressure and temperature.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        dTdP:   float
            dTdP(const. F) (K GPa-1)
        """
        _dTdPsol = self.dTdPSolidus(P)
        _dTdPliq = self.dTdPLiquidus(P)
        _F = self.F(P,T)
        if _F == 0:
            _dTdP = self.alphas/self.rhos/self.CP
        elif _F == 1:
            _dTdP = self.alphaf/self.rhof/self.CP
        else:
            _dTdP = (_F**(1/self.beta))*(_dTdPliq-_dTdPsol) + _dTdPsol
        return _dTdP

class LithologyKatz:
    """
    Lithology formatted like Katz et al. (2003) lherzolite melting model. Does
    not incorporate their hydrous melting parameterisation. Default parameter
    values taken from Katz et al. (2003) for dry lherzolite melting. All parameters
    listed below can be called.

    Parameters
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
    Mcpx:   float
        Mass fraction of cpx in the source. Controls the transition to
        low-productivity harzburgite-type melting.
    A1:     float
        Parameter used to define solidus.
    A2:     float
        Parameter used to define solidus.
    A3:     float
        Parameter used to define solidus.
    B1:     float
        Parameter used to define lherzolite-liquidus.
    B2:     float
        Parameter used to define lherzolite-liquidus.
    B3:     float
        Parameter used to define lherzolite-liquidus.
    C1:     float
        Parameter used to define liquidus.
    C2:     float
        Parameter used to define liquidus.
    C3:     float
        Parameter used to define liquidus.
    beta1:  float
        Parameter used to calculate melt fraction during cpx-present melting.
    beta2:  float
        Parameter used to calculate melt fraction during cpx-absent melting.
    r1:     float
        Parameter used to define cpx reaction coefficient.
    r2:     float
        Parameter used to define cpx reaction coefficient.
    """
    def __init__(self,Mcpx=0.15,DeltaS=300.0,A1=1085.7,A2=132.9,A3=-5.1,B1=1475.0,
                 B2=80.0,B3=-3.2,C1=1780.0,C2=45.0,C3=-2.0,beta1=1.5,beta2=1.5,
                 r1=0.5,r2=0.08,CP=1000.0,alphas=40.0,alphaf=68.0,rhos=3.3,rhof=2.9):
        self.Mcpx = Mcpx
        self.DeltaS = DeltaS
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.B1 = B1
        self.B2 = B2
        self.B3 = B3
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.beta1 = beta1
        self.beta2 = beta2
        self.r1 = r1
        self.r2 = r2
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof

    def TSolidus(self,P):
        """
        Returns the temperature of the solidus at any given pressure. Eqn(4).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tsol:   float
            Solidus temperature (degC).
        """
        _TSolidus = self.A1 + self.A2*P + self.A3*(P**2)
        return _TSolidus

    def TLiquidus(self,P):
        """
        Returns the temperature of the liquidus at any given pressure. Eqn(10).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tliq:   float
            Liquidus temperature (degC).
        """
        _TLiquidus = self.C1 + self.C2*P + self.C3*(P**2)
        return _TLiquidus

    def TLherzLiquidus(self,P):
        """
        Returns the temperature of the lherzolite liquidus at any given pressure.
        Eqn(5). This is the temperature at which the rock would be completely
        molten if cpx was remained present for the entirety of melting.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tlzliq:   float
            Lherzolite liquidus temperature (degC).
        """
        _TLherzLiquidus = self.B1 + self.B2*P + self.B3*(P**2)
        return _TLherzLiquidus

    def RescaledTcpx(self,T,P):
        """
        Calculates the rescaled temperature during cpx-present melting. (Eq 3).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (GPa).

        Returns
        -------
        T':  float
            Rescaled Temperature (dimensionless).
        """
        _TSolidus = self.TSolidus(P)
        _TLherzLiquidus = self.TLherzLiquidus(P)
        _RescaledTemperaturecpx = ((T - _TSolidus)/(_TLherzLiquidus-_TSolidus))
        return _RescaledTemperaturecpx

    def Fcpx(self,T,P):
        """
        Melt fraction during cpx-present melting at the given P and T. Eq(2).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (degC).

        Returns
        -------
        F:  float
            Melt fraction during cpx-present melting.
        """
        _RescaledT = self.RescaledTcpx(T,P)
        _Fcpx = _RescaledT**self.beta1
        return _Fcpx

    def RxnCoef(self,P):
        """
        Reaction coefficient for cpx during melting at the specified pressure. Eq(7).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        RxnCoef:    float
            Reaction Coefficient.
        """
        _RxnCoef = self.r1 + self.r2*P
        return _RxnCoef

    def FcpxOut(self,P):
        """
        Calculates the melt fraction required to exhaust cpx from the residue
        at the given pressure. Eq(6).

        Parameters
        ----------
        P:  float
            Pressure (GPa)

        Returns
        -------
        Fcpx-out:   float
            Melt fraction at which cpx is exhausted.
        """
        _FcpxOut = self.Mcpx / self.RxnCoef(P)
        return _FcpxOut

    def TcpxOut(self,P):
        """
        Calculates the temperature at which cpx will be exhausted during melting
        at the given pressure. Eq(9).

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        Tcpx-out:   float
            Temperature of cpx-exhaustion.
        """
        _TSolidus = self.TSolidus(P)
        _TcpxOut= ((self.FcpxOut(P)**(1/self.beta1)))*(self.TLherzLiquidus(P)-_TSolidus)+_TSolidus
        return _TcpxOut


    def RescaledTopx(self,T,P):
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
        T':  float
            Rescaled Temperature (dimensionless).
        """
        _TcpxOut = self.TcpxOut(P)
        _RescaledTopx = ((T-_TcpxOut)/(self.TLiquidus(P)-_TcpxOut))
        return _RescaledTopx

    def Fopx(self,T,P):
        """
        Melt fraction during cpx-absent melting at the given P and T. Eq(8).

        Parameters
        ----------
        T:  float
            Temperature (degC).
        P:  float
            Pressure (degC).

        Returns
        -------
        F:  float
            Melt fraction during cpx-absent melting.
        """
        _FcpxOut = self.FcpxOut(P)
        _FopxDry = _FcpxOut + (1-_FcpxOut)* self.RescaledTopx(T,P)**self.beta2
        return _FopxDry

    def F(self,P,T):
        """
        Wrapper for the melt fraction functions. If T and P are below the solidus,
        returns 0, if they are above the liquidus, returns 1. If below the temperature
        of cpx-exhaustion, calls the Fcpx function, otherwise calls the Fopx function.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        F:  float
            Melt fraction.
        """
        if T > self.TLiquidus(P):
            _F = 1.0
        elif T < self.TSolidus(P):
            _F = 0.0
        elif T < self.TcpxOut(P):
            _F = self.Fcpx(T,P)
        else:
            _F = self.Fopx(T,P)

        return _F

    def dTdF(self,P,T):
        """
        Calculates dT/dF(const. P). First calculates the melt fraction. If F is
        zero, returns np.inf. If F is 1, returns np.inf. Otherwise uses the
        appropriate expressions for cpx present or absent melting.

        Parameters
        ----------
        P:  float
            Pressure (GPa)
        T:  float
            Temperature (degC)

        Returns
        -------
        dTdF:   float
            dT/dF(const. P) (K).
        """
        _F = self.F(P,T)
        if _F == 0:
            dTdF = np.inf # If no melt fraction the derivative is zero. Prevents division by zero.
        elif _F < self.FcpxOut(P):
            dTdF = ((1/self.beta1))*(self.TLherzLiquidus(P)-self.TSolidus(P))*(_F**((1-self.beta1)/self.beta1))
        elif _F < 1.0:
            dTdF = ((1/self.beta2))*(self.TLiquidus(P)-self.TcpxOut(P))*(_F**((1-self.beta2)/self.beta2))
        else:
            dTdF = np.inf
        return dTdF


    def dTdP(self,P,T):
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
        dTdP:   float
            dT/dP(const. F) (K GPa-1).
        """
        _F = self.F(P,T)
        if _F == 0:
            _dTdP = self.alphas/self.rhos/self.CP
        elif _F == 1.0:
            _dTdP = self.alphaf/self.rhof/self.CP
        elif _F < self.FcpxOut(P):
            _dTdP = ((_F**(1/self.beta1))*(self.B2-self.A2+
                    2*(self.B3-self.A3)*P)+self.A2+2*self.A3*P)

        else:
            _FcpxOut = self.FcpxOut(P)
            _TLzLiq = self.TLherzLiquidus(P)
            _TSol = self.TSolidus(P)
            _dTdP1 = (((self.TLiquidus(P)-self.TcpxOut(P))/self.beta2)*
                (((_F-_FcpxOut)/(1-_FcpxOut))**((1-self.beta2)/self.beta2))*
                (((_F-_FcpxOut)/((1-_FcpxOut)**2))-(1/(1-_FcpxOut)))*
                ((self.r2*(_FcpxOut**2))/self.Mcpx))
            _dTdP2 = ((((_F-_FcpxOut)/(1-_FcpxOut))**(1/self.beta2))*
                (self.C2+2*P*self.C3-(_TLzLiq-_TSol)*
                (self.r2/(self.beta1*self.Mcpx))*(_FcpxOut**((1+self.beta1)/self.beta1))-
                (_FcpxOut**(1/self.beta1))*(self.B2+2*P*self.B3)-
                (1-(_FcpxOut**(1/self.beta1)))*(self.A2+2*P*self.A3))+
                self.A2+2*P*self.A3)
            _dTdP3 = ((self.r2/(self.beta1*self.Mcpx))*(_TLzLiq-
                _TSol)*(_FcpxOut**((1+self.beta1)/self.beta1))-
                (_FcpxOut**(1/self.beta1))*(self.B2-self.A2+2*P*(self.B3-self.A3)))
            _dTdP =  - _dTdP1 + _dTdP2 - _dTdP3

        return _dTdP


class LithologySimple:
    """
    Lithology formatted like Pertermann and Hirschmann (2002) G2 melting model.
    Default values reproduce the behaviour of G2 pyroxenite. Default thermodynamic
    constants are those used by Katz et al. (2003). All the parameters listed below
    are callable.

    Parameters
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
    a:  float
        Parameter used in calculating melt fraction.
    b:  float
        Parameter used in calculating melt fraction.
    c:  float
        Parameter used in liquidus definition.
    d:  float
        Parameter used in liquidus definition.
    e:  float
        Parameter used in solidus definition.
    f:  float
        Parameter used in solidus definition.

    """
    def __init__(self,a=0.7368,b=0.2632,c=1175,d=114,e=920,f=130,DeltaS=300.0,
                 CP=1000.0,alphas=40.0,alphaf=68.0,rhos=3.3,rhof=2.9):
        self.a=a
        self.b=b
        self.c=c
        self.d=d
        self.e=e
        self.f=f
        self.DeltaS=DeltaS
        self.CP=CP
        self.alphas=alphas
        self.alphaf=alphaf
        self.rhos=rhos
        self.rhof=rhof

    def F(self,P,T):
        """
        Calculates melt fraction at a given pressure and temperature using
        a*T'**2 + b*T', where T is the
        normalised temperature: (T-Tsolidus)/(T-Tliquidus). If P and T are
        below the the solidus, 0 is returned, if they are above the liquidus,
        1 is returned.

        Parameters
        ----------
        P:  float
            Pressure (GPa)
        T:  float
            Temperature (degC)

        Returns
        -------
        F:  float
            Melt fraction.
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        if T < _Tsol:
            _F = 0.0
        elif T > _Tliq:
            _F = 1.0
        else:
            _Tr = (T-_Tsol)/(_Tliq-_Tsol)
            _F = self.a*_Tr**2 + self.b*_Tr
        return _F

    def TLiquidus(self,P):
        """
        Calculates the liquidus temperature, at a given pressure, using
        c + d*P.

        Parameters
        ----------
        P: float
            Pressure (GPa).

        Returns
        -------
        Tliq:   float
            Liquidus temperature (degC).
        """
        _Tliq = self.c +self.d*P
        return(_Tliq)

    def TSolidus(self,P):
        """
        Calculates the solidus temperature, at a given pressure, using
        e +f*P.

        Parameters
        ----------
        P:  float
            Pressure (GPa)

        Returns
        -------
        Tsol:   float
            Solidus temperature (degC).
        """
        _Tsol = self.e +self.f*P
        return(_Tsol)

    def dTdF(self,P,T):
        """
        Calculates dT/dF(const. P) at a given pressure and temperature.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        dTdF:   float
            dT/dF(const. P) (K)
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        if T < _Tsol:
            _dTdF = np.inf
        elif T > _Tliq:
            _dTdF = np.inf
        else:
            _dTdF = (_Tliq-_Tsol)/(self.a*2*((T-_Tsol)/(_Tliq-_Tsol))+self.b)

        return _dTdF

    def dTdP(self,P,T):
        """
        Calculates dT/dP(const. F) at a given pressure and temperature.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        dTdP:   float
            dTdP(const. F) (K GPa-1)
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        _dTdP = ((self.d-self.f)*(T-_Tsol)/(_Tliq-_Tsol)+self.f)

        return _dTdP

class LithologyShorttle:
    """
    Lithology formatted like the KG1 parameterisation by Shorttle et al. (2014).
    Default values reproduce behaviour of KG1. Default thermodynamic constants
    are those used by Katz et al. (2003). All the parameters listed below are
    callable.

    Parameters
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
    A1:     float
        Constant used in solidus expression.
    A2:     float
        Constant used in solidus expression.
    A3:     float
        Constant used in solidus expression.
    B1:     float
        Constant used in cpx-out expression.
    B2:     float
        Constant used in cpx-out expression.
    B3:     float
        Constant used in cpx-out expression.
    C1:     float
        Constant used in liquidus expression.
    C2:     float
        Constant used in liquidus expression.
    C3:     float
        Constant used in liquidus expression.
    a:  float
        Constant used in cpx-present melt fraction expression.
    b:  float
        Constant used in cpx-present melt fraction expression.
    c:  float
        Constant used in cpx-absent melt fraction expression.
    d:  float
        Constant used in cpx-absent melt fraction expression.
    alpha:  float
        Exponent used in the cpx-present melt fraction expression.
    beta:   float
        Exponent used in the cpx-absent melt fraction expression.
    """
    def __init__(self,DeltaS=300.0,CP=1000.0,alphas=40.0,alphaf=68.0,rhos=3.3,rhof=2.9,
                 A1=1095.4,A2=124.1,A3=-4.7,B1=1179.6,B2=157.2,B3=-11.1,C1=1780.0,C2=45.0,C3=-2.0,
                 a=0.3187864,b=0.4154,c=0.7341864,d=0.2658136,alpha=2,beta=1.5):
         self.DeltaS = DeltaS
         self.CP = CP
         self.alphas = alphas
         self.alphaf = alphaf
         self.rhos = rhos
         self.rhof = rhof
         self.A1 = A1
         self.A2 = A2
         self.A3 = A3
         self.B1 = B1
         self.B2 = B2
         self.B3 = B3
         self.C1 = C1
         self.C2 = C2
         self.C3 = C3
         self.a = a
         self.b = b
         self.c = c
         self.d = d
         self.alpha = alpha
         self.beta = beta

    def TSolidus(self,P):
        """
        Returns solidus temperature at a given pressure. T = A1 + A2*P + A3*P**2.

        Parameters
        ----------
        P:  float, or list of floats, or array of floats
            Pressure in GPa.

        Returns
        -------
        T:  float, or list of floats, or array of floats
            Solidus Temperature in degC.
        """
        _T = self.A1 + self.A2*P + self.A3*P**2
        return _T

    def TLiquidus(self,P):
        """
        Returns liquidus temperature at a given pressure. T = C1 + C2*P + C3*P**2.

        Parameters
        ----------
        P:  float, or list of floats, or array of floats
            Pressure in GPa.

        Returns
        -------
        T: float, or list of floats, or array of floats
            Liquidus temperature in degC.
        """
        _T = self.C1 + self.C2*P + self.C3*P**2
        return _T

    def TCpxOut(self,P):
        """
        Returns the temperature of cpx-exhaustion at a given pressure. T = B1 + B2*P + B3*P**2.

        Parameters
        ----------
        P: float, or list of floats, or array of floats
            Pressure in GPa.

        Returns
        -------
        T: float, or list of floats, or array of floats
            Cpx-exhaustion temperature in degC.
        """
        _T = self.B1 + self.B2*P + self.B3*P**2
        return _T

    def dTdP(self,P,T):
        """
        Returns dT/dP (constant F) at a given pressure and temperature.

        Parameters
        ----------
        P:  float
            Pressure in GPa.
        T:  float
            Temperature in degC.

        Returns
        -------
        dTdP:   float
            dT/dP (constant F) in K GPa-1.
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        _Tcpx = self.TCpxOut(P)

        if T < _Tcpx:
            _dTdP = -(-(T-_Tsol)/(_Tcpx-_Tsol)*(self.B2+2*self.B3*P-self.A2-self.A3*2*P) -
                    self.A2 - self.A3*2*P)
        else:
            _dTdP = -(-(T-_Tcpx)/(_Tliq-_Tcpx)*(self.C2+self.C3*2*P-self.B2-self.B3*2*P) -
                    self.B2 - 2*self.B3*P)

        return _dTdP

    def dTdF(self,P,T):
        """
        Returns dT/dF (constant P) at a given pressure and temperature. If below
        the solidus, or above the liquidus, np.inf is returned.

        Parameters
        ----------
        P:  float
            Pressure in GPa.
        T:  float
            Temperature in degC.

        Returns
        -------
        dTdF:   float
            dT/dF (constant P) in K.
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        _Tcpx = self.TCpxOut(P)

        if T < _Tsol:
            _dTdF = np.inf
        elif T < _Tcpx:
            _dTdF = (_Tcpx-_Tsol)/(self.b + self.a*self.alpha*((T-_Tsol)/(_Tcpx-_Tsol))**(self.alpha-1))
        elif T < _Tliq:
            _dTdF = (_Tliq - _Tcpx)/(self.d*self.beta*((T-_Tcpx)/(_Tliq-_Tcpx))**(self.beta-1))
        else:
            _dTdF = np.inf

        return _dTdF

    def F(self,P,T):
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
        P:  float
            Pressure in GPa.
        T:  float
            Temperature in degC.

        Returns
        -------
        F:  float
            Melt fraction between 0 and 1.
        """
        _Tsol = self.TSolidus(P)
        _Tliq = self.TLiquidus(P)
        _Tcpx = self.TCpxOut(P)

        if T < _Tsol:
            _F = 0.0
        elif T > _Tliq:
            _F = 1.0
        elif T < _Tcpx:
            Tf = (T-_Tsol)/(_Tcpx-_Tsol)
            _F = self.a*Tf**self.alpha + self.b*Tf
        else:
            Tf = (T-_Tcpx)/(_Tliq-_Tcpx)
            _F = self.d*Tf**self.beta + self.c
        return _F

class LithologyNonMelting:
    """
    Material that does not melt, i.e. Harzburgite in Shorttle et al. (2014) and
    Matthews et al. (2016). Default thermodynamic constants are those used by
    Katz et al. (2003).

    Properties
    ----------
    DeltaS:
        Entropy change during melting. Set to 0.0.
    CP:
        Heat capacity (J Kg-1 K-1)
    alphas:
        Thermal expansivity of the solid (K-1)
    alphaf:
        Thermal expansivity of the melt. Set to 0.0.
    rhos:
        Density of the solid (g cm-3).
    rhof:
        Density of the melt. Set to 0.0.

    Parameters
    ----------
    CP:     float
        Heat capacity (J Kg-1 K-1)
    alphas:     float
        Thermal expansivity of the solid (K-1)
    rhos:   float
        Density of the solid (g cm-3)
    """
    def __init__(self,CP=1000.0,alphas=40.0,rhos=3.3):
         self.DeltaS = 0.0
         self.CP = CP
         self.alphas = alphas
         self.alphaf = 0.0
         self.rhos = rhos
         self.rhof = 0.0

    def F(self,P,T):
        """
        Melt Fraction. Returns 0.0.

        Parameters
        ----------
        P:
            Pressure. There to maintain consistancy within lithology definitions.
        T:
            Temperature. There to maintain consistancy within lithology definitions.
        """
        return 0.0

    def dTdF(self,P,T):
        """
        dTdF(constP). Returns np.inf.

        Parameters
        ----------
        P:
            Pressure. There to maintain consistancy within lithology definitions.
        T:
            Temperature. There to maintain consistancy within lithology definitions.
        """
        return np.inf

    def dTdP(self,P,T):
        """
        dTdP(constF). Returns 0.0.

        Parameters
        ----------
        P:
            Pressure. There to maintain consistancy within lithology definitions.
        T:
            Temperature. There to maintain consistancy within lithology definitions.
        """
        return 0.0

    def TSolidus(self,P):
        """
        Solidus temperature. Returns np.inf.

        Parameters
        ----------
        P:
            Pressure. There to maintain consistancy within lithology definitions.
        """
        return np.inf

    def TLiquidus(self,P):
        """
        Liquidus temperature. Returns np.inf

        Parameters
        ----------
        P:
            Pressure. There to maintain consistancy within lithology definitions.
        """
        return np.inf

class mantle:
    """
    Mantle class consists of one or more lithology classes, in a particular proportion.
    Properties that change with pressure and temperature are never made a part of the
    class. Callable properties are:

    number_lithologies: int
        the number of lithologies in the mantle class
    CP: list of floats
        the heat capacities of the lithologies
    alphaf: list of floats
        the thermal expansion coefficients of the melts produced by each lithology.
    alphas: list of floats
        the thermal expansion coefficients of each lithology
    rhof: list of floats
        the densities of the melts produced by each lithology
    rhos: list of floats
        the densities of each lithology
    DeltaS: list of floats
        the entropy change on melting of each lithology.


    Parameters
    ----------
    lithologies:    list of lithology objects
        A list of the defined lithology objects.
    proportions:    list of floats
        The mass ratios of the lithologies, doesn't need to be normalised.
    names:  list of strings
        The names of the lithologies. If False, default names will be chosen.


    """
    def __init__(self,lithologies,proportions,names=False):
        self.lithologies = lithologies
        if isinstance(proportions,list):
            proportions = np.array(proportions)
        self.proportions = proportions/sum(proportions)
        self.number_lithologies = np.shape(self.lithologies)[0]
        if names==False:
            names = list()
            for i in range(self.number_lithologies):
                names.append('default '+str(i))
        self.names = names

        self.CP = np.zeros(self.number_lithologies)
        self.alphaf = np.zeros(self.number_lithologies)
        self.alphas = np.zeros(self.number_lithologies)
        self.rhof = np.zeros(self.number_lithologies)
        self.rhos = np.zeros(self.number_lithologies)
        self.DeltaS = np.zeros(self.number_lithologies)
        for i in range(self.number_lithologies):
            self.CP[i] = self.lithologies[i].CP
            self.alphaf[i] = self.lithologies[i].alphaf
            self.alphas[i] = self.lithologies[i].alphas
            self.rhof[i] = self.lithologies[i].rhof
            self.rhos[i] = self.lithologies[i].rhos
            self.DeltaS[i] = self.lithologies[i].DeltaS

    def bulk_properties(self,P=False,T=False):
        """
            Calculates the bulk thermodynamic properties of the solid or partially
            molten mantle.

            Parameters
            ----------
            P:  float or bool
                The pressure of interest. If False the properties of the solid mantle
                will be returned instead.
            T:  float or bool
                The temperature of interest. If False the properties of the solid mantle
                will be returned instead.

            Returns
            -------
            bulk_constants:     pandas Series
                The bulk alpha, CP and rho for the mantle at the given P and T, labelled
                as such. Uses pandas Series rather than a dict in order to allow access via
                self.bulk_properties().alpha etc.
        """
        _F = np.zeros(self.number_lithologies)
        if T != False:
            for i in range(self.number_lithologies):
                _F[i] = self.lithologies[i].F(P,T)
        _alpha = sum(self.proportions*self.alphaf*_F+self.proportions*self.alphas*(1-_F))
        _CP = sum(self.proportions*self.CP)
        _rho = sum(self.proportions*self.rhof*_F+self.proportions*self.rhos*(1-_F))

        return pd.Series({'alpha':_alpha,'CP':_CP,'rho':_rho})

    def solidus_intersection(self,Tp):
        """
        Finds the pressure at which each lithology's solidus will be intersected,
        assuming the mantle follows the solid adiabat up until that point.

        Parameters
        ----------
        Tp:     float
            The mantle potential temperature in degC.

        Returns
        -------
        intersection:   np.array
            The pressure of solidus intersection of each lithology (in the order
            passed when creating the mantle class).
        """
        _intersect = np.zeros(self.number_lithologies)
        for i in range(self.number_lithologies):
            def f_solve(P):
                return self.lithologies[i].TSolidus(P) - self.adiabat(P,Tp)
            _intersect[i] = fsolve(f_solve,3.0)[0]
        return _intersect


    def solidus_intersection_isobaric(self,P):
        """
        Finds the pressure at which each lithology's solidus will be intersected,
        assuming the mantle is heated isobarically.

        Parameters
        ----------
        P:     float
            The pressure of interest in GPa

        Returns
        -------
        intersection:   np.array
            The temperature of solidus intersection of each lithology (in the order
            passed when creating the mantle class).
        """
        _intersect = np.zeros(self.number_lithologies)
        for i in range(self.number_lithologies):
            _intersect[i] = self.lithologies[i].TSolidus(P)
        return _intersect


    def adiabat(self,P,Tp):
        """
        Calculates the temperature of the solid mantle at a given pressure, knowing
        the potential temperature.

        Parameters
        ----------
        P:  float or list or np.array
            Pressure in GPa.
        Tp:     float
            Potential temperature in degC.

        Returns
        -------
        T:  float or list or np.array
            Temperature of the mantle at the given pressure and Tp.
        """
        _T = ((Tp+273.15)*np.exp(self.bulk_properties(P).alpha/
              (self.bulk_properties(P).rho*self.bulk_properties(P).CP)*P) - 273.15)
        return _T

    def F(self,P,T):
        """
        Calculates the melt fraction of each lithology at a given pressure and
        temperature. Acts as a lookup function.

        Parameters
        ----------
        P:  float
            Pressure in GPa
        T:  float
            Temperature in degC

        Returns
        -------
        F:  np.array
            Array containing the melt fraction of each lithology, in the order
            the lithologies were passed to the mantle class.
        """
        _F = np.zeros(self.number_lithologies)
        for i in range(self.number_lithologies):
            _F[i] = self.lithologies[i].F(P,T)
        return _F

    def dFdP(self,P,T):
        """
        Calculates the value of dFdP for each lithology at the given pressure
        and temperature, using Eq(26) of Phipps Morgan (2001).

        Parameters
        ----------
        P:  float
            Pressure in GPa.
        T:  float
            Temperature in degC.

        Returns
        -------
        dFdP    np.array
            Array of dFdP values for each lithology, in the order in which they
            were passed to the mantle class.
        """
        _dFdP = np.zeros(self.number_lithologies)

        _dTdP = np.zeros(self.number_lithologies)
        _dTdF = np.zeros(self.number_lithologies)
        _F = np.zeros(self.number_lithologies)

        for i in range(self.number_lithologies):
            _dTdP[i] = self.lithologies[i].dTdP(P,T)
            _dTdF[i] = self.lithologies[i].dTdF(P,T)
            _F[i] = self.lithologies[i].F(P,T)


        _lithologies_melting = np.where((_F>0)&(_F<1))[0]

        if np.shape(_lithologies_melting)[0] > 0:
            for i in range(np.shape(_lithologies_melting)[0]):
                _not_key = [True]*self.number_lithologies
                _key = _lithologies_melting[i]
                _not_key[_key] = False


                # Equation (26) from PM2001 to find dFdP of first lithology
                _top = (self.bulk_properties(P,T).CP/(T+273.15)*_dTdP[_key] -
                        self.bulk_properties(P,T).alpha/self.bulk_properties(P,T).rho +
                        sum(self.proportions[_not_key]*self.DeltaS[_not_key]*
                            (_dTdP[_key]-_dTdP[_not_key])/_dTdF[_not_key]))

                _bottom = (self.proportions[_key]*self.DeltaS[_key] +
                            sum(self.proportions[_not_key]*self.DeltaS[_not_key]*_dTdF[_key]/_dTdF[_not_key]) +
                            self.bulk_properties(P,T).CP/(T+273.15)*_dTdF[_key])

                _dFdP[_key] = - _top/_bottom

        return _dFdP

    def adiabatic_gradient(self,P,T):
        """
        Calculates dTdP if melting has gone to completion (or hasn't started) for
        the bulk mantle.

        Parameters
        ----------
        P:  float
            Pressure in GPa.
        T:  float
            Temperature in degC.

        Returns
        -------
        dTdP:   float
                The adiabatic gradient in C/GPa
        """

        _bulk = self.bulk_properties(P,T)

        _dTdP = _bulk['alpha']*(T+273.15)/_bulk['rho']/_bulk['CP']

        return _dTdP


    def dTdP(self,P,T,dFdP):
        """
        Calculates dTdP using Eq(28) of Phipps Morgan (2001). Picks the lithology
        to use by the one with the largest increase in melt fraction with decompression
        (though this choice shouldn't matter). This is something that should be used to
        check the behaviour of the code.

        Parameters
        ----------
        P:  float
            Pressure in GPa. Needed to find dTdP(const F), and dTdF(const P).
        T:  float
            Temperature in degC. Needed to find dTdP(const F), and dTdF(const P).
        dFdP:   np.array
            The return from the dFdP function. Array of the dFdP values of each
            lithology and the same pressure and temperature as passed to the
            function. This function doesn't recall the dFdP function in order
            to save re-calculating the same numbers multiple times.

        Returns
        -------
        dTdP:   float
            The thermal gradient in the melting region at the P, and T of interest.
        """
        _melting_lithologies = np.where(dFdP<0)[0]
        if np.shape(_melting_lithologies)[0] > 0:
            _key = np.argmin(dFdP[dFdP<0])
            _key = _melting_lithologies[_key]
            _dTdP = self.lithologies[_key].dTdP(P,T) + self.lithologies[_key].dTdF(P,T)*dFdP[_key]
        else:
            _dTdP = self.lithologies[0].dTdP(P,T)

        return _dTdP

    def IsobaricMelt_1D(self,Tstart,P,dT=0.1):
        """
        Calculates the amount of melt generated, and the mantle's temperature, after
        an interval of melting occuring due to mantle being instantaneously placed
        above its solidus.

        The intention of this function was to handle high Tp cases where the solidus
        is always exceeded, not to produce physically meaningful results, but to
        allow the tails of Tp distributions when inverting to be reasonably approximated.

        Parameters
        ----------
        Tstart:     float
            The temperature (degC) at which to place the solid mantle.
        P:  float
            The pressure at which to perform the calculation.
        dT:     float
            The interval of discretisation of temperature increase from the solidus.
        """

        _solidus_intersection = self.solidus_intersection_isobaric(P)

        # Calculate the entropy lost associated with cooling solid material to the
        # solidus temperature
        _solT = np.nanmin(_solidus_intersection)
        _DeltaS_cool = - self.bulk_properties()['CP']*np.log((_solT+273)/(Tstart+273))

        _DeltaS_melt = 0
        _T = _solT + dT
        while _DeltaS_melt < _DeltaS_cool and _T < Tstart:
            _DeltaS_melt = (np.sum(self.F(P,_T)*self.proportions*self.DeltaS) +
                            - self.bulk_properties(P,_T)['CP']*np.log((_solT+273)/(_T+273)))
            _T = _T + dT

        return _T




    def AdiabaticMelt_1D(self,Tp,Pstart=8.0,Pend=0.01,steps=1001,ReportSSS=True):
        """
        Performs simultaneous integration of dFdP and dTdP in order to obtain
        the thermal gradient through the melting region. F of each lithology is then
        calculated using the P,T path. Integration is performed using a 4th order
        Runge-Kutta algorithm.

        The T-P path is allowed to overstep the solidus on the step prior to the
        start of melting.

        Parameters
        ----------
        Tp:     float
            The potential temperature (degC) at which to perform the calculation.
        Pstart:     float
            The pressure (in GPa) at which to begin upwelling. Default is 8 GPa.
        Pend:   float
            The pressure (in GPa) at which to stop upwelling. Default is 0 GPa.
        steps:  int
            The number of dP increments to split the melting region into. Default
            is 1001.
        ReportSSS:  bool
            Print to the console if the start is above the solidus of one of the
            lithologies. Either way the code will calculate the melt fraction at
            this point by conserving entropy. Set to False if this is deliberate.

        Returns
        -------
        MeltingColumn:  MeltingColumn_1D object
            The results are returned in a 1D Melting Column object, further
            calculations, e.g. crustal thickness may then be performed on this
            object, if desired.
        """
        _T = np.zeros(steps)
        _T[0] = self.adiabat(Pstart,Tp)
        _P = np.linspace(Pstart,Pend,steps)
        _dP = (Pend-Pstart)/(steps-1)
        _F = np.zeros([steps,self.number_lithologies])

        if _T[0] > np.nanmin(self.solidus_intersection_isobaric(Pstart)):
            if ReportSSS == True:
                print('WARNING! SUPER SOLIDUS START')
            _T[0] = self.IsobaricMelt_1D(_T[0],Pstart)

        for i in range(steps):
            if i == 0:
                _F[i] = self.F(_P[0],_T[0])
            else:
                if np.shape(np.where((_F[i-1]>0))[0])[0] == 0:
                    _T[i] = self.adiabat(_P[i],Tp)
                    _F[i] = self.F(_P[i],_T[i])
                elif np.shape(np.where((_F[i-1]>0)&(_F[i-1]<1))[0])[0] == 0:
                    _j1 = self.adiabatic_gradient(_P[i-1],_T[i-1])
                    _j2 = self.adiabatic_gradient(_P[i-1]+_dP/2,_T[i-1]+_dP/2*_j1)
                    _j3 = self.adiabatic_gradient(_P[i-1]+_dP/2,_T[i-1]+_dP/2*_j2)
                    _j4 = self.adiabatic_gradient(_P[i],_T[i-1]+_dP*_j3)

                    _T[i] = _T[i-1] + _dP/6*(_j1+2*_j2+2*_j3+_j4)
                    _F[i] = self.F(_P[i],_T[i])
                else:
                    _k1 = self.dFdP(_P[i-1],_T[i-1])
                    _j1 = self.dTdP(_P[i-1],_T[i-1],_k1)
                    _k2 = self.dFdP(_P[i-1]+_dP/2,_T[i-1]+_dP/2*_j1)
                    _j2 = self.dTdP(_P[i-1]+_dP/2,_T[i-1]+_dP/2*_j1,_k2)
                    _k3 = self.dFdP(_P[i-1]+_dP/2,_T[i-1]+_dP/2*_j2)
                    _j3 = self.dTdP(_P[i-1]+_dP/2,_T[i-1]+_dP/2*_j2,_k3)
                    _k4 = self.dFdP(_P[i],_T[i-1]+_dP*_j3)
                    _j4 = self.dTdP(_P[i],_T[i-1]+_dP*_j3,_k4)

                    _T[i] = _T[i-1] + _dP/6*(_j1+2*_j2+2*_j3+_j4)
                    _F[i] = self.F(_P[i],_T[i])

        _results = pd.DataFrame(_F,columns=self.names)
        _results['P'] = _P
        _results['Temperature'] = _T

        return MeltingColumn_1D(_results,self,Tp)

    def PlotBoundaries(self,Pmax=8.0,Pmin=0.0,steps=1000,T_F=1600,show=True):
        """
        Generates 2 plots, one showing the P-T relationship of the solidii and
        liquidii of each of the lithologies, and one which shows the melt fraction
        of each lithology at fixed temperature as a function of pressure.

        Parameters
        ----------
        Pmax:   float
            The maximum pressure (GPa) to display. Default is 8 GPa.
        Pmin:   float
            The minimum pressure (GPa) to display. Default is 0 GPa.
        steps:  float
            The discretization to use. Default is 1000.
        T_F:    float
            The temperature (degC) at which to calculate melt fractions.
        show:   bool
            Display the plot, or not.

        Returns
        -------
        f:  matplotlib.figure
            The figure object, for adjustments or saving.
        """
        _P = np.linspace(Pmin,Pmax,steps)

        f,a = plt.subplots(1,2,sharey='row')

        for i in range(self.number_lithologies):
            _Tsol = self.lithologies[i].TSolidus(_P)
            _Tliq = self.lithologies[i].TLiquidus(_P)
            _F = list()
            for j in range(steps):
                _F.append(self.lithologies[i].F(_P[j],T_F))

            if isinstance(self.lithologies[i],LithologyNonMelting) == False:
                a[0].plot(_Tsol,_P,label=self.names[i]+' solidus')
                a[0].plot(_Tliq,_P,label=self.names[i]+' liquidus')
                a[1].plot(_F,_P,label=self.names[i]+' at '+str(T_F)+' C')

            if (isinstance(self.lithologies[i],LithologyKatz) or
                isinstance(self.lithologies[i],LithologyKLB1) or
                isinstance(self.lithologies[i],LithologyKG1)):
                _TLzLiq = self.lithologies[i].TLherzLiquidus(_P)
                a[0].plot(_TLzLiq,_P,label=self.names[i]+' lz liquidus')
            if isinstance(self.lithologies[i],LithologyShorttle):
                _Tcpx = self.lithologies[i].TCpxOut(_P)
                a[0].plot(_Tcpx,_P,label=self.names[i]+' cpx out')


        a[0].legend()
        a[1].legend()
        a[0].set_xlabel('T ($^\circ$C)')
        a[1].set_xlabel('F')
        a[0].set_ylabel('P (GPa)')
        a[0].invert_yaxis()

        a[0].tick_params('x',labeltop=True,labelbottom=False)
        a[0].xaxis.set_label_position('top')
        a[1].tick_params('x',labeltop=True,labelbottom=False)
        a[1].xaxis.set_label_position('top')

        if show == True:
            plt.show()
        return f

class MeltingColumn_1D():
    """
    Class for storing the results of a 1D multi-lithology melting model.

    Parameters
    ----------
    calculation_results:    pandas.DataFrame
        Dataframe with columns 'P' for Pressure in GPa, 'Temperature' for Temperature in
        degrees C, Remaining columns for melt fraction from each lithology.
    mantle:     mantle object
        The mantle object used to generate the melting column.
    Tp:     float
        The potential temperature used to generate the melting column, if applicable.
    """
    def __init__(self,calculation_results,mantle,Tp=False):
        self.P = calculation_results.P
        self.Temperature = calculation_results.Temperature
        _cols = calculation_results.columns.tolist()
        _cols.remove('P')
        _cols.remove('Temperature')
        self.F = calculation_results[_cols]
        self.mantle = mantle
        self.Tp = Tp
        self.F_total = np.zeros(np.shape(self.F)[0])
        for i in range(self.mantle.number_lithologies):
            self.F_total = self.F_total + self.mantle.proportions[i]*self.F[self.mantle.names[i]]

    def plot(self,solidii=True,show=True):
        """
        Generates a standard plot showing the thermal gradient and melt fractions of
        each lithology.

        Parameters
        ----------
        solidii     bool
            Plot solidii, or not.
        show    bool
            Show the plot after generating it, or not.

        Returns
        -------
        f   matplotlib.figure
            The generated figure.
        """
        f,a = plt.subplots(1,2,sharey='row')
        _lith = self.F.columns
        for i in range(np.shape(_lith)[0]):
            a[1].plot(self.F.iloc[:,i],self.P,label=_lith[i])
        a[0].plot(self.Temperature,self.P,label='Thermal Gradient',c='k')
        a[1].plot(self.F_total,self.P,label='total',c='k',ls='--')
        if solidii == True and self.mantle!=False:
            _P = np.linspace(np.min(self.P),np.max(self.P),1000)
            for i in range(self.mantle.number_lithologies):
                if isinstance(self.mantle.lithologies[i],LithologyNonMelting) == False:
                    _T = self.mantle.lithologies[i].TSolidus(_P)
                    a[0].plot(_T,_P,label=self.mantle.names[i]+' solidus')

        if self.Tp != False:
            a[0].text(0.95,0.95,'T$_p$ = '+str(self.Tp)+' $^\circ$C',
                     transform=a[0].transAxes,va='top',ha='right')

        a[0].invert_yaxis()

        a[0].set_ylabel('Pressure (GPa)')
        a[0].set_xlabel('Temperature ($^\circ$C)')
        a[1].set_xlabel('Melt Fraction')

        a[0].legend(loc='lower left')
        a[1].legend()

        a[0].tick_params('x',labeltop=True,labelbottom=False)
        a[0].xaxis.set_label_position('top')
        a[1].tick_params('x',labeltop=True,labelbottom=False)
        a[1].xaxis.set_label_position('top')

        if show == True:
            plt.show()

        return f

    def integrate_tri(self):
        """
        Perform an integration over the melting region, assuming it is triangular and
        passively upwelling. Adds the following attributes to the class:

        dtcdP:  array of floats
            the results from Eq6 for the total melt fraction.
        tc_int:     array of floats
            integrated crustal thickness as a function of pressure (up to 0 GPa)
        tc_P_int:   array of floats
            the pressure exerted by the integrated crustal thickness as a function
            of pressure (up to 0 GPa).
        tc:     float
            The integrated crustal thickness at the point where the pressure it
            exerts is equal to the calculation pressure.
        P_base_of_crust:    float
            The pressure at the base of the crust, at the point where the pressure
            the generated crust exerts is equal to the calculation pressure.
        tc_lithology_contributions_int:     pandas DataFrame
            The integrated proportion of generated crust derived from each lithology
            as a function of pressure.
        tc_lithology_contributions:     pandas Series
            The integrated proportion of generated crust derived from each lithology
            at the pressure where P(calculation) = P(exerted by generated crust)

        Returns
        -------
        tc:     float
            The crustal thickness where the pressure exerted from the produced crust
            equals the pressure of melting.

        """
        _rho = self.mantle.bulk_properties().rho
        _g=9.81
        _tc = 1/(_rho*_g*1e3)*self.F_total/(1-self.F_total)
        _tc_lith = 1/(_rho*_g*1e3)*self.F*self.mantle.proportions/(1-np.tile(self.F_total,[np.shape(self.F)[1],1]).T)
        _tc_int = np.zeros(np.shape(self.P)[0])
        _tc_lith_int = np.zeros(np.shape(_tc_lith))
        _tc_intP = np.zeros(np.shape(self.P)[0])
        _tc_found = False
        _P_basecrust = False
        _tc_lith_found = False
        for i in range(np.shape(self.P)[0]):
            if i != 0:
                _tc_int[i] = _tc_int[i-1]+ _tc[i]*np.abs(self.P[i]-self.P[i-1])
                _tc_lith_int[i] = _tc_lith_int[i-1] + _tc_lith.iloc[i]*np.abs(self.P[i]-self.P[i-1])
                _tc_intP[i] = _tc_int[i]*_rho*_g*1e3
                if _tc_intP[i] > self.P[i] and _tc_found == False:
                    _tc_found = _tc_int[i]
                    _P_basecrust = self.P[i]
                    _tc_lith_found = _tc_lith_int[i]


        self.dtcdP = _tc
        self.tc_int = _tc_int*1e6
        self.tc_P_int = _tc_intP
        self.tc = _tc_found*1e6
        self.P_base_of_crust = _P_basecrust
        self.tc_lithology_contributions_int = _tc_lith/_tc_lith.sum()
        self.tc_lithology_contributions = _tc_lith_found/sum(_tc_lith_found)

        return _tc_found*1e6

    def MeltCrystallisationT(self,ShallowMeltP=False,MeltStorageP=False,liqdTdP=39.16):
        """
        Identifies the crystallisation temperature of the deepest and shallowest melts,
        according to the technique used by Matthews et al. (2016).

        ShallowMeltP:   float
            The pressure (in GPa) at which the shallowest melt should be extracted. If set
            to False (as is default) this will be taken as the base of the crust. If triangular
            integration has not been done, this will result in an error.
        MeltStorageP:   float
            The pressure at which crystallisation is happening. If set to False (as is default),
            the base of the crust will be used. If triangular integration has not been done,
            this will result in an error.
        liqdTdP:    float
            The clapeyron slope of the liquidus (K GPa-1), the default value is 39.16,
            from equation (15) of Putirka  (2008).
        """


        # Set crystallisation Pressure
        if MeltStorageP == False:
            MeltStorageP = self.P_base_of_crust
        if ShallowMeltP == False:
            ShallowMeltP = self.P_base_of_crust

        self.T_crystallisation = self.Temperature - (self.P-MeltStorageP)*liqdTdP

        # Extract crystallisation temperatures of interest
#        FirstMeltRow = r.F_total[r.F_total>0].idxmin()
        self.DeepMeltTcrys = {}
        for l in self.mantle.names:
            if (self.F[l]>0).any():
                FirstMeltRow = self.F[l][self.F[l]>0].idxmin()
                self.DeepMeltTcrys[l] = self.T_crystallisation.iloc[FirstMeltRow]
            else:
                self.DeepMeltTcrys[l] = np.nan
        self.TcrysMax = np.nanmax(list(self.DeepMeltTcrys.values()))

        LastMeltRow = self.P[self.P>ShallowMeltP].idxmin()
        self.TcrysMin = self.T_crystallisation.iloc[LastMeltRow]

        return self.TcrysMin, self.TcrysMax
