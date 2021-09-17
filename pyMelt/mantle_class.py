import numpy as np
from warnings import warn
import pandas as pd
from scipy.optimize import fsolve

from pyMelt.meltingcolumn_classes import MeltingColumn_1D


class mantle:
    """
    The mantle class consists of one or more lithology classes, in a particular proportion. The
    mantle class contains the methods used for doing melting calculations and for calculating the
    properties of a heterogeneous mantle.

    Parameters
    ----------
    lithologies :    list of lithology objects
        A list of lithology instances.
    proportions :    list or numpy.array of floats
        The mass ratios of the lithologies, doesn't need to be normalised.
    names :  list of strings or None, default: None
        The names of the lithologies. If False, default names will be chosen.

    Attributes
    ----------
    number_lithologies : int
        the number of lithologies in the mantle class
    CP :     list of floats
        the heat capacities of the lithologies
    alphaf : list of floats
        the thermal expansion coefficients of the melts produced by each lithology.
    alphas : list of floats
        the thermal expansion coefficients of each lithology
    rhof :   list of floats
        the densities of the melts produced by each lithology
    rhos :   list of floats
        the densities of each lithology
    DeltaS : list of floats
        the entropy change on melting of each lithology.

    """
    def __init__(self, lithologies, proportions, names=None):
        self.lithologies = lithologies
        if isinstance(proportions, list):
            proportions = np.array(proportions)
        self.proportions = proportions/sum(proportions)
        self.number_lithologies = np.shape(self.lithologies)[0]
        if names is None:
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

    def bulk_properties(self, P=None, T=None):
        """
            Calculates the bulk thermodynamic properties of the solid or partially
            molten mantle.

            Parameters
            ----------
            P : float or None, default: None
                The pressure of interest. If None, the properties of the solid mantle will be
                returned.
            T : float or None, default: None
                The temperature of interest. If None, the properties of the solid mantle will be
                returned.

            Returns
            -------
            dict
                The bulk alpha, CP and rho for the mantle at the given P and T, labelled
                as such.
        """
        F = np.zeros(self.number_lithologies)
        if T is not None:
            for i in range(self.number_lithologies):
                F[i] = self.lithologies[i].F(P, T)
        alpha = sum(self.proportions*self.alphaf*F + self.proportions*self.alphas*(1-F))
        CP = sum(self.proportions*self.CP)
        rho = sum(self.proportions*self.rhof*F + self.proportions*self.rhos*(1-F))

        return {'alpha': alpha, 'CP': CP, 'rho': rho}

    def solidus_intersection(self, Tp):
        """
        Finds the pressure at which each lithology's solidus will be intersected,
        assuming the mantle follows the solid adiabat up until that point.

        Parameters
        ----------
        Tp : float
            The mantle potential temperature in degC.

        Returns
        -------
        numpy.array
            The pressure of solidus intersection of each lithology.
        """
        intersect = np.zeros(self.number_lithologies)
        for i in range(self.number_lithologies):
            def f_solve(P):
                return self.lithologies[i].TSolidus(P) - self.adiabat(P, Tp)
            # Check that the lithology actually melts:
            if self.lithologies[i].TSolidus(3.0) is np.inf:
                intersect[i] = np.nan
            else:
                try:
                    intersect[i] = fsolve(f_solve, 3.0)[0]
                except Exception:
                    intersect[i] = np.nan
        return intersect

    def solidus_intersection_isobaric(self, P):
        """
        Finds the pressure at which each lithology's solidus will be intersected,
        assuming the mantle is heated isobarically.

        Parameters
        ----------
        P : loat
            The pressure of interest in GPa

        Returns
        -------
        np.array
            The temperature of solidus intersection of each lithology.
        """
        intersect = np.zeros(self.number_lithologies)
        for i in range(self.number_lithologies):
            intersect[i] = self.lithologies[i].TSolidus(P)
        return intersect

    def adiabat(self, P, Tp):
        """
        Calculates the actual temperature of the solid mantle at a given pressure, given the
        potential temperature.

        Parameters
        ----------
        P :  float or numpy.array
            Pressure in GPa.
        Tp : float or numpy.array
            Potential temperature in degC.

        Returns
        -------
        float or numpy.array
            Temperature of the mantle at the given pressure and Tp.
        """
        bulk_props = self.bulk_properties(P)
        T = ((Tp + 273.15) * np.exp(bulk_props['alpha'] /
             (bulk_props['rho'] * bulk_props['CP']) * P) - 273.15)
        return T

    def F(self, P, T):
        """
        Calculates the melt fraction of each lithology at a given pressure and
        temperature.

        Parameters
        ----------
        P : float
            Pressure in GPa
        T : float
            Temperature in degC

        Returns
        -------
        np.array
            Array containing the melt fraction of each lithology.
        """
        F = np.zeros(self.number_lithologies)
        for i in range(self.number_lithologies):
            F[i] = self.lithologies[i].F(P, T)
        return F

    def dFdP(self, P, T):
        """
        Calculates the value of dFdP for each lithology at the given pressure
        and temperature, using Eq(26) of Phipps Morgan (2001).

        Parameters
        ----------
        P : float
            Pressure in GPa.
        T : float
            Temperature in degC.

        Returns
        -------
        np.array
            Array of dFdP values for each lithology.
        """
        dFdP = np.zeros(self.number_lithologies)

        dTdP = np.zeros(self.number_lithologies)
        dTdF = np.zeros(self.number_lithologies)
        F = np.zeros(self.number_lithologies)

        for i in range(self.number_lithologies):
            dTdP[i] = self.lithologies[i].dTdP(P, T)
            dTdF[i] = self.lithologies[i].dTdF(P, T)
            F[i] = self.lithologies[i].F(P, T)

        lithologies_melting = np.where((F > 0) & (F < 1))[0]

        if np.shape(lithologies_melting)[0] > 0:
            for i in range(np.shape(lithologies_melting)[0]):
                not_key = [True]*self.number_lithologies
                key = lithologies_melting[i]
                not_key[key] = False

                bulk_props = self.bulk_properties(P, T)

                # Equation (26) from PM2001 to find dFdP of first lithology
                top = (bulk_props['CP'] / (T+273.15) * dTdP[key] -
                       bulk_props['alpha']/bulk_props['rho'] +
                       sum(self.proportions[not_key] * self.DeltaS[not_key] *
                           (dTdP[key] - dTdP[not_key]) / dTdF[not_key]))

                bottom = (self.proportions[key] * self.DeltaS[key] +
                          sum(self.proportions[not_key] * self.DeltaS[not_key] *
                              dTdF[key] / dTdF[not_key]) +
                          bulk_props['CP'] / (T+273.15) * dTdF[key])

                dFdP[key] = - top/bottom

        return dFdP

    def adiabatic_gradient(self, P, T):
        """
        Calculates dTdP if melting has gone to completion (or hasn't started) for the bulk mantle.

        Parameters
        ----------
        P : float
            Pressure in GPa.
        T : float
            Temperature in degC.

        Returns
        -------
        float
            The adiabatic gradient in C/GPa
        """

        bulk_props = self.bulk_properties(P, T)

        dTdP = bulk_props['alpha'] * (T+273.15) / bulk_props['rho'] / bulk_props['CP']

        return dTdP

    def dTdP(self, P, T, dFdP=None):
        """
        Calculates dTdP using Eq(28) of Phipps Morgan (2001). Picks the lithology to use by the
        one with the largest increase in melt fraction with decompression (though this choice
        shouldn't matter).

        Parameters
        ----------
        P :    float
            Pressure in GPa.
        T :    float
            Temperature in degC.
        dFdP : np.array or None, default: None
            The value of dFdP at the same T and P. If None, the dFdP method will be called. In
            most melting calculations the value will already have been calculated, so passing the
            value will save computational time.

        Returns
        -------
        float
            The thermal gradient in the melting region at the P and T of interest.
        """
        if dFdP is None:
            dFdP = self.dFdP(P, T)

        melting_lithologies = np.where(dFdP < 0)[0]

        if np.shape(melting_lithologies)[0] > 0:
            key = np.argmin(dFdP[dFdP < 0])
            key = melting_lithologies[key]
            dTdP = self.lithologies[key].dTdP(P, T) + self.lithologies[key].dTdF(P, T) * dFdP[key]
        else:
            dTdP = self.adiabatic_gradient(P, T)

        return dTdP

    def IsobaricMelt_1D(self, Tstart, P, dT=0.1):
        """
        Calculates the amount of melt generated, and the mantle's temperature, after
        an interval of melting occuring due to mantle being instantaneously placed
        above its solidus.

        The intention of this function was to handle high Tp cases where the solidus
        is always exceeded, not to produce physically meaningful results, but to
        allow the tails of Tp distributions to be approximated reasonably when inverting.

        Parameters
        ----------
        Tstart : float
            The temperature (degC) at which to place the solid mantle.
        P :      float
            The pressure at which to perform the calculation.
        dT :     float
            The interval of discretisation of temperature increase from the solidus.

        Returns
        -------
        float
            Temperature of mantle following the melting step.
        """

        solidus_intersection = self.solidus_intersection_isobaric(P)

        # Calculate the entropy lost associated with cooling solid material to the
        # solidus temperature
        solT = np.nanmin(solidus_intersection)
        bulk_props = self.bulk_properties(P, solT)
        DeltaS_cool = - bulk_props['CP'] * np.log((solT + 273.15) / (Tstart + 273.15))

        DeltaS_melt = 0
        T = solT + dT
        while DeltaS_melt < DeltaS_cool and T < Tstart:
            bulk_props = self.bulk_properties(P, T)
            DeltaS_melt = (np.sum(self.F(P, T) * self.proportions * self.DeltaS) +
                           - bulk_props['CP'] * np.log((solT + 273.15) / (T + 273.15)))
            T = T + dT

        return T

    def AdiabaticMelt_1D(self, Tp, Pstart=8.0, Pend=0.01, steps=101, ReportSSS=True):
        """
        Performs simultaneous integration of dFdP and dTdP to obtain the thermal gradient
        through the melting region. F of each lithology is then calculated using the P,T path.
        Integration is performed using a 4th order Runge-Kutta algorithm.

        The T-P path is allowed to overstep the solidus on the step prior to the start of melting.

        Parameters
        ----------
        Tp :        float
            The potential temperature (degC) at which to perform the calculation.
        Pstart :    float, default: 8.0
            The pressure (in GPa) at which to begin upwelling.
        Pend :      float, default: 0.0
            The pressure (in GPa) at which to stop upwelling.
        steps :     int, default: 101
            The number of dP increments to split the melting region into.
        ReportSSS : bool, default: True
            Print to the console if the start is above the solidus of one of the lithologies.
            Either way the code will calculate the melt fraction at this point by conserving
            entropy. Set to False if you don't want to be warned.

        Returns
        -------
        pyMelt.MeltingColumn_1D
            The results are returned in a 1D Melting Column instance, further calculations, e.g.
            crustal thickness may then be performed on this instance, if desired.
        """

        T = np.zeros(steps)
        T[0] = self.adiabat(Pstart, Tp)
        P = np.linspace(Pstart, Pend, steps)
        dP = (Pend - Pstart) / (steps - 1)
        F = np.zeros([steps, self.number_lithologies])

        if T[0] > np.nanmin(self.solidus_intersection_isobaric(Pstart)):
            if ReportSSS is True:
                warn("Super solidus start")
            T[0] = self.IsobaricMelt_1D(T[0], Pstart)

        for i in range(steps):
            if i == 0:
                F[i] = self.F(P[0], T[0])
            else:
                # If melting has not started
                if np.shape(np.where((F[i-1] > 0))[0])[0] == 0:
                    T[i] = self.adiabat(P[i], Tp)
                    F[i] = self.F(P[i], T[i])

                # If melting has gone to completion
                elif np.shape(np.where((F[i-1] > 0) & (F[i-1] < 1))[0])[0] == 0:
                    j1 = self.adiabatic_gradient(P[i-1], T[i-1])
                    j2 = self.adiabatic_gradient(P[i-1] + dP/2, T[i-1] + dP/2*j1)
                    j3 = self.adiabatic_gradient(P[i-1] + dP/2, T[i-1] + dP/2*j2)
                    j4 = self.adiabatic_gradient(P[i], T[i-1] + dP*j3)

                    T[i] = T[i-1] + dP/6*(j1 + 2*j2 + 2*j3 + j4)
                    F[i] = self.F(P[i], T[i])

                # If melting is ongoing
                else:
                    k1 = self.dFdP(P[i-1], T[i-1])
                    j1 = self.dTdP(P[i-1], T[i-1], k1)
                    k2 = self.dFdP(P[i-1] + dP/2, T[i-1] + dP/2*j1)
                    j2 = self.dTdP(P[i-1] + dP/2, T[i-1] + dP/2*j1, k2)
                    k3 = self.dFdP(P[i-1] + dP/2, T[i-1] + dP/2*j2)
                    j3 = self.dTdP(P[i-1] + dP/2, T[i-1] + dP/2*j2, k3)
                    k4 = self.dFdP(P[i], T[i-1] + dP*j3)
                    j4 = self.dTdP(P[i], T[i-1] + dP*j3, k4)

                    T[i] = T[i-1] + dP/6 * (j1 + 2*j2 + 2*j3 + j4)
                    F[i] = self.F(P[i], T[i])

        results = pd.DataFrame(F, columns=self.names)
        results['P'] = P
        results['T'] = T

        return MeltingColumn_1D(results, self, Tp)
