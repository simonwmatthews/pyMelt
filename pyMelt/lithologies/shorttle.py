from pyMelt.lithology_class import lithology

import numpy as np

class kg1(lithology):
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
