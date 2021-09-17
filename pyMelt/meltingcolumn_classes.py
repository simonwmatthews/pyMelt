"""
===============
Melting Columns
===============

the meltingcolumn_classes module provides the melting column classes. At the moment it consists of
a single melting column- a simple 1D melting column.

"""

import numpy as np
from pyMelt.core import InputError
import matplotlib.pyplot as plt


class MeltingColumn_1D():
    """
    Class for storing the results of a 1D multi-lithology melting model.

    Parameters
    ----------
    calculation_results : pandas.DataFrame
        Dataframe with columns 'P' for Pressure in GPa, 'T' for Temperature in
        degrees C, Remaining columns for melt fraction from each lithology.
    mantle : pyMelt.Mantle
        The mantle class used to generate the melting column.
    Tp : float
        The potential temperature used to generate the melting column, if applicable.

    Attributes
    ----------
    P : pandas.Series
        The pressures of the melting column in GPa.
    T : pandas.Series
        The temperatures of the melting column in degC at each pressure.
    F : pandas.DataFrame
        The melt fractions of each lithology at each pressure.
    mantle : pyMelt.Mantle
        The mantle class used to generate the melting model
    Tp : float
        The potential temperature (in degC) used to generate the melts.
    F_total : numpy.array
        The total melt fraction at each pressure.
    dtcdP : numpy.array or None
        The results from Eq6 for the total melt fraction following integration. Available following
        integration.
    tc_int : numpy.array or None
        Integrated crustal thickness as a function of pressure (up to 0 GPa). Available following
        integration.
    tc_P_int : numpy.array or None
        The pressure exerted by the integrated crustal thickness as a function of pressure (up to
        0 GPa). Available following integration.
    tc : float or None
        The integrated crustal thickness at the point where the pressure it exerts is equal to the
        calculation pressure. Available following integration.
    P_base_of_crust : float or None
        The pressure at the base of the crust, at the point where the pressure the generated crust
        exerts is equal to the calculation pressure. Available following integration.
    tc_lithology_contributions_int : pandas.DataFrame or None
        The integrated proportion of generated crust derived from each lithology as a function of
        pressure. Available following integration.
    tc_lithology_contributions : pandas.Series or None
        The integrated proportion of generated crust derived from each lithology at the pressure
        where P(calculation) = P(exerted by generated crust). Available following integration.
    """

    def __init__(self, calculation_results, mantle, Tp=None):
        self.P = calculation_results.P
        self.T = calculation_results['T']

        cols = calculation_results.columns.tolist()
        cols.remove('P')
        cols.remove('T')
        self.F = calculation_results[cols]

        self.mantle = mantle
        self.Tp = Tp
        self.F_total = np.zeros(np.shape(self.F)[0])

        for i in range(self.mantle.number_lithologies):
            self.F_total = self.F_total + self.mantle.proportions[i]*self.F[self.mantle.names[i]]

        self.dtcdP = None
        self.tc_int = None
        self.tc_P_int = None
        self.tc = None
        self.P_base_of_crust = None
        self.tc_lithology_contributions_int = None
        self.tc_lithology_contributions = None

    def plot(self):
        """
        Generates a plot showing the thermal gradient and melt fractions ofeach lithology.

        Returns
        -------
        (matplotlib.figure, matplotlib.axes)
            The generated figure and axes objects.
        """
        f, a = plt.subplots(1, 2, sharey='row', dpi=100)

        lith = self.F.columns

        for i in range(np.shape(lith)[0]):
            a[1].plot(self.F.iloc[:, i], self.P, label=lith[i])

        a[0].plot(self.T, self.P, label='Thermal Gradient', c='k')
        a[1].plot(self.F_total, self.P, label='Total', c='k', ls='--')

        P = np.linspace(np.min(self.P), np.max(self.P), 101)
        for i in range(self.mantle.number_lithologies):
            T = self.mantle.lithologies[i].TSolidus(P)
            a[0].plot(T, P, label=self.mantle.names[i]+' solidus')

        if self.Tp is not None:
            a[0].text(0.95, 0.95, 'T$_p$ = {:.0f} °C'.format(self.Tp),
                      transform=a[0].transAxes, va='top', ha='right')

        a[0].invert_yaxis()

        a[0].set_ylabel('Pressure (GPa)')
        a[0].set_xlabel('Temperature (°C)')
        a[1].set_xlabel('Melt Fraction')

        a[0].legend(loc='lower left')
        a[1].legend()

        a[0].tick_params('x', labeltop=True, labelbottom=False)
        a[0].xaxis.set_label_position('top')
        a[1].tick_params('x', labeltop=True, labelbottom=False)
        a[1].xaxis.set_label_position('top')

        return f, a

    def integrate_tri(self, P_base_existingLith=0.0, extract_melt=False):
        """
        Perform an integration over the melting region, assuming it is triangular and passively
        upwelling.

        Parameters
        ----------
        P_base_existingLith : float, default: 0.0
            The pressure at the base of any pre-existing lithosphere. The calculated thickness will
            be added to this value. Set to 0.0 for mid-ocean ridges, set to non-zero for
            continental rifting.

        extract_melt : bool, default: False
            If set to True, the melts will be extracted from the system. If False (default) the
            melt produced is added to the top of the melting column, and the calculation will stop
            once the pressure exerted from the newly made crust is equal to the pressure in the
            calculation step.

        Returns
        -------
        float
            The crustal thickness where the pressure exerted from the produced crust equals the
            pressure of melting.

        """
        rho = self.mantle.bulk_properties()['rho']
        g = 9.81
        tc = 1.0 / (rho * g * 1e3) * self.F_total / (1.0 - self.F_total)
        tc_lith = (1.0 / (rho * g * 1e3) * self.F * self.mantle.proportions
                   / (1.0 - np.tile(self.F_total, [np.shape(self.F)[1], 1]).T))
        tc_int = np.zeros(np.shape(self.P)[0])
        tc_lith_int = np.zeros(np.shape(tc_lith))
        tc_intP = np.zeros(np.shape(self.P)[0])
        tc_found = False
        P_basecrust = False
        tc_lith_found = False

        for i in range(np.shape(self.P)[0]):
            if i != 0:
                tc_int[i] = tc_int[i-1] + tc[i] * np.abs(self.P[i] - self.P[i-1])
                tc_lith_int[i] = (tc_lith_int[i-1]
                                  + tc_lith.iloc[i] * np.abs(self.P[i] - self.P[i-1]))
                tc_intP[i] = tc_int[i] * rho * g * 1e3
                if(extract_melt is False and tc_intP[i] + P_base_existingLith > self.P[i]
                   and tc_found is False):
                    tc_found = tc_int[i]
                    P_basecrust = self.P[i]
                    tc_lith_found = tc_lith_int[i]
                elif (extract_melt is True
                      and (i == np.shape(self.P)[0] - 1 or P_base_existingLith > self.P[i])
                      and tc_found is False):
                    tc_found = tc_int[i]
                    P_basecrust = self.P[i]
                    tc_lith_found = tc_lith_int[i]

        self.dtcdP = tc
        self.tc_int = tc_int * 1e6
        self.tc_P_int = tc_intP
        self.tc = tc_found * 1e6
        self.P_base_of_crust = P_basecrust
        self.tc_lithology_contributions_int = tc_lith / tc_lith.sum()
        self.tc_lithology_contributions = tc_lith_found / sum(tc_lith_found)

        return tc_found * 1e6

    def MeltCrystallisationT(self, ShallowMeltP=None, MeltStorageP=None, liqdTdP=39.16):
        """
        Identifies the crystallisation temperature of the deepest and shallowest melts,
        according to the technique used by Matthews et al. (2016).

        Parameters
        ----------
        ShallowMeltP : float or None, default: None
            The pressure (in GPa) at which the shallowest melt should be extracted. If set
            to None (as is default) this will be taken as the base of the crust. If integration
            has not been performed, this will result in an error.
        MeltStorageP : float, default:  None
            The pressure at which crystallisation is happening. If set to False (as is default),
            the base of the crust will be used. If triangular integration has not been done,
            this will result in an error.
        liqdTdP : float, default: 39.16
            The clapeyron slope of the liquidus (K GPa-1), the default value is 39.16,
            from equation (15) of Putirka  (2008).

        Returns
        -------
        float
            The minimum crystallisation temperature (degC)
        float
            The maximum crystallisation temperature (degC)
        """

        # Set crystallisation Pressure
        if MeltStorageP is None:
            if self.P_base_of_crust is None:
                raise InputError("You must perform an integration before calculating Tcrys.")
            MeltStorageP = self.P_base_of_crust
        if ShallowMeltP is None:
            ShallowMeltP = self.P_base_of_crust

        self.T_crystallisation = self.T - (self.P - MeltStorageP) * liqdTdP

        self.DeepMeltTcrys = {}
        for l in self.mantle.names:
            if (self.F[l] > 0).any():
                FirstMeltRow = self.F[l][self.F[l] > 0].idxmin()
                self.DeepMeltTcrys[l] = self.T_crystallisation.iloc[FirstMeltRow]
            else:
                self.DeepMeltTcrys[l] = np.nan
        self.TcrysMax = np.nanmax(list(self.DeepMeltTcrys.values()))

        LastMeltRow = self.P[self.P > ShallowMeltP].idxmin()
        self.TcrysMin = self.T_crystallisation.iloc[LastMeltRow]

        return self.TcrysMin, self.TcrysMax
