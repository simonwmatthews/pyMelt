#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')


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

    def integrate_tri(self, P_base_existingLith=0.0, extract_melt = False):
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

        Parameters
        ----------
        P_base_existingLith:    float
            The pressure at the base of any pre-existing lithosphere. The calculated
            thickness will be added to this value. Set to 0.0 for mid-ocean ridges,
            set to non-zero for continental rifting.

        extract_melt:       bool
            If set to True, the melts will be extracted from the system. If False (default)
            the melt produced is added to the top of the melting column, and the calculation
            will stop once the pressure exerted from the newly made crust is equal to the
            pressure in the calculation step.


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
                if extract_melt == False and _tc_intP[i] + P_base_existingLith > self.P[i] and _tc_found == False:
                    _tc_found = _tc_int[i]
                    _P_basecrust = self.P[i]
                    _tc_lith_found = _tc_lith_int[i]
                elif extract_melt == True and (i == np.shape(self.P)[0] - 1 or P_base_existingLith > self.P[i]) and _tc_found == False:
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
