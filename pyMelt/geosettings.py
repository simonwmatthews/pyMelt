"""
===================
Geological Settings
===================

The geological settings module provides the classes required for extracting setting-specific
information from the melting columns, e.g., the crustal thickness and average melt composition
at a spreading centre.

"""

import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
from copy import copy
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from warnings import warn

import pyMelt.chemistry

# Weighting function
def weighting_expdecay(P, weighting_wavelength, weighting_amplitude=1.0):
    """
    Weights melts according to an exponential decay function, with the greatest weighting at the
    base of the melting region, and zero weighting at the top.

    Parameters
    ----------
    P : numpy.array or pandas.Series
        The complete set of pressure values for the melting region.
    weighting_wavelength : float
        The wavelength of the exponential decay (applied to pressures normalised to vary between
        0 and 1).
    weighting_amplitude : float, default: 1.0
        The amplitude of the exponential function.
    """
    Pmax = _np.max(MeltingColumn.P)
    Pmax = _np.min(MeltingColumn.P)

    w = amplitude * _np.exp(- wavelength * (Pmax - P) / (Pmax - Pmin))

    return w

class GeoSetting(object):
    """
    Base clase for geological settings.

    Parameters
    ----------
    MeltingColumn : pyMelt.meltingcolumn_classes.MeltingColumn
        The melting column from which to construct the geological setting.
    weightingFunction : function or None, default: None
        A function used to apply an additional weighting to melts during homogenisation, perhaps
        for simulating the behaviour of active upwelling. The result from the function will be
        added to the triangular weighting applied already. The function must take the melting
        region pressures (as a numpy array or pandas Series) as its first argument, any other
        arguments will be passed in kwargs.


    Attributes
    ----------
    MeltingColumn : pyMelt.meltingcolumn_classes.MeltingColumn
        The melting column from which the geological setting was constructed
    lithologies : dict
        Dictionary containing the DataFrames of the states of each lithology during melting.
    mantle : pyMelt.Mantle
        The mantle object from which the melting column was calculated.
    """

    def __init(self, MeltingColumn, weightingFunction = None, **kwargs):
        self.MeltingColumn = MeltingColumn
        self.lithologies = MeltingColumn.lithologies
        self.mantle = MeltingColumn.mantle
        self.weightingFunction = weightingFunction
        self.kwargs = kwargs

    def crystallisation_chemistry(self, mineralProportions, fractionate=True,
                                  D = pyMelt.chemistry.defaultD):
        """
        Calculates the concentrations of elements in the homogenised magma following an interval
        of crystallisation.

        Parameters
        ----------
        fractionate : bool, default: True
            If True, the calculation will follow fractional crystallisation.
        mineralProportions: dict
            The proportions of each mineral, relative to the system total (1.0). The minerals must
            have the same name as used in the D argument. If the default partition coefficients
            are used the mineral labels are olv, cpx, opx, plg, spn, grt.
        D : pandas.DataFrame, default: pyMelt.chemistry.defaultD
            The partition coefficients for each element and mineral. Columns are minerals, rows
            are elements.

        Returns
        -------
        pandas.Series
            The concentrations of the trace elements in the evolved melt.
        """

        X = 1 - sum(mineralProportions.values())
        Dbulk = _pd.Series(index = self.chemistry.index, data = [0]*len(self.chemistry))
        for mineral in mineralProportions:
            Dbulk += mineralProportions[mineral] * D[mineral] / (1 - X)
        if fractionate is True:
            evolved_chemistry = self.chemistry * X ** (Dbulk - 1.0)
        else:
            evolved_chemistry = self.chemistry / (X * (1 - Dbulk) + Dbulk)
        return evolved_chemistry

    def plot_spider(self, normalisation = 'PM', plot_instantaneous=False, plot_original=True,
                    crystal_fraction=None, element_order=None, **kwargs):
        """
        Plot a basic spider diagram of the chemical composition, optionally alongside the
        instantaneous melts and/or an evolved melt.

        Parameters
        ----------
        normalisation : dict or str, default: 'PM'
            How to normalise the values. To use the built in default options pass a string:
            - 'PM' gives the Primitive Mantle of Palme and O'Neill (2013)
            - 'CI' gives the chondritic composition of Palme and O'Neill (2013)
            - 'DM' gives the depleted mantle composition of Workman and Hart (2005)

            Otherwise, pass a dict of elements (keys) and their concentrations.
        plot_instantaneous : bool, default: False
            Plot the range of the instantaneous melts (or accumulated melts at each step) as a
            field.
        plot_original : bool, default: True
            Plot the original homogenised chemistry.
        crystal_fraction : None or dict, default: None
            Plot an evolved melt having crystallised crystals in the proportions passed. If using
            non-default options for GeoSetting.crystallisation_chemistry() pass these as additional
            arguments
        element_order : None or list, default: None
            Use to adjust the order of the elements, otherwise they will be given in the order of
            the chemistry attribute.

        Returns
        -------
        matplotlib.figure, matplotlib.axes

        """
        f, a = _plt.subplots(dpi=150)

        default_norms = {'PM': pyMelt.chemistry.palme13_pm,
                         'CI': pyMelt.chemistry.palme13_ci,
                         'DM': pyMelt.chemistry.workman05_dmm}

        if isinstance(normalisation, str):
            normalisation = default_norms[normalisation]

        if element_order is None:
            element_order = list(self.chemistry.keys())

        if plot_instantaneous is True:
            normed_hi = []
            normed_lo = []
            for el in element_order:
                hi = 0.0
                lo = 1e16
                for lith in self.MeltingColumn.lithologies:
                    lithmax = _np.max(self.MeltingColumn.lithologies[lith][el])
                    lithmin = _np.min(self.MeltingColumn.lithologies[lith][el])
                    if lithmax > hi:
                        hi = lithmax
                    if lithmin < lo:
                        lo = lithmin
                    normed_hi.append(hi/normalisation[el])
                    normed_lo.append(lo/normalisation[el])
            a.fill_between(range(len(element_order)), normed_lo, normed_hi, alpha=0.2)

        if plot_original is True:
            normed = []
            for el in element_order:
                normed.append(self.chemistry[el]/normalisation[el])
            a.plot(range(len(element_order)), normed)

        if crystal_fraction is not None:
            unnormed = self.crystallisation_chemistry(crystal_fraction, **kwargs)
            normed = []
            for el in element_order:
                normed.append(unnormed[el]/normalisation[el])
            a.plot(range(len(element_order)), normed)

        a.set_yscale('log')
        a.set_xticks(range(len(element_order)))
        a.set_xticklabels(element_order)
        a.set_ylabel('Normalised concentration')

        return f, a

    def _weighting_coefficients(self, P):
        """
        Calculates the weighting coefficients, or returns 0 if no weighting function is supplied.

        Parameters
        ----------
        P : numpy.array or pandas.Series
            The complete set of pressure values for the melting region.
        """
        if self.weightingFunction is not None:
            weights = self.weightingFunction(P, self.kwargs)
        else:
            weights = _np.zeros(_np.shape(P))

        return weights


class SpreadingCentre(GeoSetting):
    """
    Implementation of a spreading centre, representing either mid-ocean ridge spreading or
    continental rift spreading. The crustal thickness is calculated assuming a triangular melting
    region.

    Parameters
    ----------
    MeltingColumn : pyMelt.meltingcolumn_classes.MeltingColumn
        The melting column from which to construct the geological setting.
    P_lithosphere : float, default: 0.0
        The pressure at the base of the lithosphere in a continental rift. If this includes the
        igneous crust, set `extract_melt` to True. Defaults to 0.0, the case of a a mid-ocean
        ridge.
    extract_melt : bool, default: False
        Should integration continue beyond the calculated igneous crustal thickness (+
        lithosphere)? Set to False for modelling a mid-ocean ridge.
    steps: int or None, default: 10001
        The number of steps over which to perform the crustal thickness calculation. If None, the
        same number of steps as the melting calculation will be used. Otherwise, the melting
        results will be interpolated over the new number of steps. This primarily controls how
        finely the crustal thickness is determined, i.e., the size of the P overstep in obtaining
        the crustal thickness.
    weightingFunction : function or None, default: None
        A function used to apply an additional weighting to melts during homogenisation, perhaps
        for simulating the behaviour of active upwelling. The result from the function will be
        added to the triangular weighting applied already. The function must take the melting
        region pressures (as a numpy array or pandas Series) as its first argument, any other
        arguments will be passed in kwargs.

    Attributes
    ----------
    MeltingColumn : pyMelt.meltingcolumn_classes.MeltingColumn
        The melting column from which the geological setting was constructed
    lithologies : dict
        Dictionary containing the DataFrames of the states of each lithology during melting.
    mantle : pyMelt.Mantle
        The mantle object from which the melting column was calculated.
    P_lithosphere : float
        The pressure at the base of the lithosphere.
    tc : float
        The crustal thickness, in km.
    P_base_of_crust : float
        The pressure at the base of the crust, in GPa.
    lithology_contributions : pandas.Series
        The relative contributions of each lithology to the pooled melt.
    chemistry : pandas.Series
        The homogenised melt composition
    """

    def __init__(self, MeltingColumn, P_lithosphere=0.0, extract_melt=False, steps=10001,
                 weightingFunction=None, **kwargs):
        self.MeltingColumn = MeltingColumn
        self.lithologies = copy(MeltingColumn.lithologies)
        self.mantle = MeltingColumn.mantle
        self.P = copy(self.MeltingColumn.P)
        self.T = copy(self.MeltingColumn.T)
        self.F = copy(self.MeltingColumn.F)
        self.P_lithosphere = P_lithosphere
        self.weightingFunction = weightingFunction
        self.kwargs = kwargs

        # These will be set by the `integrate_tri` method.
        self.tc = None
        self.P_base_of_crust = None
        self.lithology_contributions = None

        # This will be set by the _homogenise_chemistry method.
        self.chemistry = None

        # Do the integrations for a spreading centre
        self._integrate_tri(P_lithosphere, extract_melt, steps=steps)
        self._homogenise_chemistry()

        # Remove melts from the column that would have never been produced
        if extract_melt is False:
            self.F = self.F[self.P > self.P_base_of_crust]
            self.T = self.T[self.P > self.P_base_of_crust]
            for lith in self.lithologies:
                self.lithologies[lith] = self.lithologies[lith][self.P > self.P_base_of_crust]
            self.P = self.P[self.P > self.P_base_of_crust]
        else:
            self.F = self.F[self.P > self.P_lithosphere]
            self.T = self.F[self.P > self.P_lithosphere]
            for lith in self.lithologies:
                self.lithologies[lith] = self.lithologies[lith][self.P > self.P_lithosphere]
            self.P = self.F[self.P > self.P_lithosphere]

    def _integrate_tri(self, P_base_existingLith=0.0, extract_melt=False, steps=1001):
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

        # Calculate the contributions to tc (without additional weighting)
        tc = 1.0 / (rho * g * 1e3) * self.MeltingColumn.F / (1.0 - self.MeltingColumn.F)

        # Do the same for the individual lithologies
        Flith = _pd.DataFrame()
        for i in range(self.mantle.number_lithologies):
            Flith[self.mantle.names[i]] = self.lithologies[self.mantle.names[i]].F

        tc_lith = (1.0 / (rho * g * 1e3) * Flith * self.mantle.proportions
                   / (1.0 - _np.tile(self.F, [self.mantle.number_lithologies, 1]).T))

        if steps is None:
            P = self.P
            tc_lith = tc_lith.to_numpy()
        else:
            interp_f = interp1d(self.P, tc)
            interp_lith = []
            for i in range(self.mantle.number_lithologies):
                interp_lith.append(interp1d(self.P, tc_lith[self.mantle.names[i]]))
            P = _np.linspace(_np.nanmax(self.P), _np.nanmin(self.P), steps)
            tc = _np.zeros(_np.shape(P)[0])
            for i in range(len(P)):
                tc[i] = interp_f(P[i])
            tc_lith = _np.zeros([_np.shape(P)[0], self.mantle.number_lithologies])
            for i in range(self.mantle.number_lithologies):
                for j in range(_np.shape(P)[0]):
                    tc_lith[j, i] = interp_lith[i](P[j])

        weights = self._weighting_coefficients(P)

        tc_int = _np.zeros(_np.shape(P)[0])
        tc_lith_int = _np.zeros([_np.shape(P)[0], _np.shape(tc_lith)[1]])
        tc_intP = _np.zeros(_np.shape(P)[0])
        tc_found = False
        P_basecrust = False
        tc_lith_found = False

        for i in range(_np.shape(P)[0]):
            if i != 0:
                tc_int[i] = tc_int[i - 1] + tc[i] * (_np.abs(P[i] - P[i - 1]) + weights)
                tc_lith_int[i] = (tc_lith_int[i - 1]
                                  + tc_lith[i] * (_np.abs(P[i] - P[i - 1]) + weights))
                tc_intP[i] = tc_int[i] * rho * g * 1e3
                if(extract_melt is False and tc_intP[i] + P_base_existingLith > P[i]
                   and tc_found is False):
                    tc_found = tc_int[i]
                    P_basecrust = P[i]
                    tc_lith_found = tc_lith_int[i]
                elif (extract_melt is True
                      and (i == _np.shape(P)[0] - 1 or P_base_existingLith > P[i])
                      and tc_found is False):
                    tc_found = tc_int[i]
                    P_basecrust = P[i]
                    tc_lith_found = tc_lith_int[i]

        self.tc = tc_found * 1e6
        self.P_base_of_crust = P_basecrust
        lith_cont = tc_lith_found / sum(tc_lith_found)
        self.lithology_contributions = {}
        for i, lith in zip(range(self.mantle.number_lithologies), self.mantle.names):
            self.lithology_contributions[lith] = lith_cont[i]

    def _homogenise_chemistry(self):
        """
        Homogenises the melt compositions according to a triangular melting region.
        """
        # Find the elements which all melting lithologies have:
        first_lithology = True
        species = []
        for lith in self.mantle.names:
            if first_lithology is False and self.lithology_contributions[lith] > 0:
                specieslith = list(self.lithologies[lith].columns)[3:]
                species = [s for s in species if s in specieslith]
            elif first_lithology is True and self.lithology_contributions[lith] > 0:
                species = list(self.lithologies[lith].columns)[3:]
                first_lithology = False

        weights = self._weighting_coefficients(self.MeltingColumn.P)

        cm = _np.zeros([len(species)])
        for lith in self.mantle.names:
            if self.lithology_contributions[lith] > 0:
                # Get the chemistry for the lithology
                c = self.lithologies[lith][species]
                # Convert it a numpy array
                c = c.to_numpy()
                # Change nans to 0.0 to avoid affecting summation
                c = _np.nan_to_num(c, nan=0.0)
                # Get the melt fractions for the lithology
                f = self.lithologies[lith].F.to_numpy()
                # If instantaneous melts, need to pool over columns first
                for j in range(len(species)):
                    if self.MeltingColumn._species_calc_type[lith][j] == 'instantaneous':
                        cnormed = _np.zeros(_np.shape(c)[0])
                        df = f[1:] - f[: - 1]
                        for i in range(_np.shape(c)[0] - 1):
                            if f[i] > 0 and f[i - 1] > 0:
                                cnormed[i + 1] = (_np.sum(0.5 * (c[1:i + 1, j] + c[0:i, j])
                                                     * df[:i], axis=0) / f[i])
                            elif f[i] > 0:
                                cnormed[i + 1] = (_np.sum(c[1:i + 1, j]
                                                     * df[:i], axis=0) / f[i])
                            else:
                                cnormed[i + 1] = 0
                        c[:, j] = cnormed
                # Normalise melts for this lithology
                c = trapz((1 + weights) * c * f[:, None] / (1.0 - f[:, None]),
                          self.lithologies[lith]['Pressure'].to_numpy()[:], axis=0)
                c = c / trapz((1 + weights) * f / (1 - f),
                              self.lithologies[lith]['Pressure'].to_numpy())

                # Normalise melts for all lithologies
                cm += c * self.lithology_contributions[lith]
            self.chemistry = _pd.Series(cm, species)

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
            MeltStorageP = self.P_base_of_crust
        if ShallowMeltP is None:
            ShallowMeltP = self.P_base_of_crust

        T_crystallisation = self.T - (self.P - MeltStorageP) * liqdTdP

        DeepMeltTcrys = {}
        for l in self.mantle.names:
            if (self.lithologies[l].F > 0).any():
                FirstMeltRow = self.lithologies[l].F[self.lithologies[l].F > 0].idxmin()
                DeepMeltTcrys[l] = T_crystallisation.iloc[FirstMeltRow]
            else:
                DeepMeltTcrys[l] = _np.nan
        TcrysMax = _np.nanmax(list(DeepMeltTcrys.values()))

        LastMeltRow = self.P[self.P > ShallowMeltP].idxmin()
        TcrysMin = T_crystallisation.iloc[LastMeltRow]

        return TcrysMin, TcrysMax


class OceanIsland(GeoSetting):
    """
    Implementation of an ocean island, representing mantle upwelling beneath lithosphere. The
    melt flux is calculated assuming flow in a deformable plume conduit (Turcotte and Schubert,
    2002). At present a constant rate of decompression throughout the conduit is assumed, likely
    leading to inaccuracies in the estimated melt chemistry.

    Parameters
    ----------
    MeltingColumn : pyMelt.meltingcolumn_classes.MeltingColumn
        The melting column from which to construct the geological setting.
    P_lithosphere : float, default: 0.0
        The pressure at the base of the lithosphere in a continental rift. If this includes the
        igneous crust, set `extract_melt` to True. Defaults to 0.0, the case of a a mid-ocean
        ridge.
    relative_density : float or None, default: None
        The value of (ambient-density - plume-density) in kg m-3.
    viscosity : float, default: 1e19
        The viscosity of the mantle plume in Pa s. Default value is 1e19 Pa s, after Shorttle et
        al. (2014).
    radius : float, default: 1e5
        The plume radius in m. Default is 1e5 m (or 100 km), after Shorttle et al. (2014).
    weighting_function : function or None
        A function of pressure allowing non-uniform weighting of melts throughout the melting
        region. Useful for simulating active upwelling, for example.

    Attributes
    ----------
    MeltingColumn : pyMelt.meltingcolumn_classes.MeltingColumn
        The melting column from which the geological setting was constructed
    lithologies : dict
        Dictionary containing the DataFrames of the states of each lithology during melting.
    mantle : pyMelt.Mantle
        The mantle object from which the melting column was calculated.
    P_lithosphere : float
        The pressure at the base of the lithosphere.
    melt_flux : float or None
        If the melt flux has been calculated, it will be stored here, in m3 s-1.
    lithology_contributions : pandas.Series
        The relative contributions of each lithology to the pooled melt.
    chemistry : pandas.Series
        The homogenised melt composition
    """

    def __init__(self, MeltingColumn, P_lithosphere, relative_density=None, viscosity=1e19,
                 radius=1e5):
        self.MeltingColumn = MeltingColumn
        self.lithologies = copy(MeltingColumn.lithologies)
        self.P_lithosphere = P_lithosphere
        self.mantle = MeltingColumn.mantle
        self.P = copy(self.MeltingColumn.P)
        self.T = copy(self.MeltingColumn.T)
        self.F = copy(self.MeltingColumn.F)

        # Remove melts that are produced more shallow than the base of the lithosphere
        self.F = self.F[self.P > self.P_lithosphere]
        self.T = self.T[self.P > self.P_lithosphere]
        for lith in self.lithologies:
            self.lithologies[lith] = self.lithologies[lith][self.P > self.P_lithosphere]
        self.P = self.P[self.P > self.P_lithosphere]

        # Extract the lithology contributions:
        self.lithology_contributions = _pd.Series({})
        for lith in self.mantle.names:
            self.lithology_contributions[lith] = (_np.nanmax(self.lithologies[lith].F)
                                                  / _np.nanmax(self.F))

        # Calculate the melt flux.
        if relative_density is not None:
            self.melt_flux = self.calc_melt_flux(relative_density, viscosity, radius)
        else:
            self.melt_flux = None

        self.chemistry = None
        self._homogenise_chemistry()

    def calc_melt_flux(self, relative_density, viscosity, radius):
        r"""
        Calculates the melt flux for given plume conduit parameters. Assumes a deformable conduit
        with constant upwelling velocity, after Turcotte & Schubert (2002).

        Parameters
        ----------
        relative_density : float
            The value of (plume-density - ambient-density) in kg m-3.
        viscosity : float
            The viscosity of the mantle plume in Pa s. Default value is 1e19 Pa s, after Shorttle
            et al. (2014).
        radius : float
            The plume radius in m. Default is 1e5 m (or 100 km), after Shorttle et al. (2014).

        Notes
        -----
        The volume flux is obtained from the equation:

        .. math::

           Q_v = \frac{\pi}{8} \frac{\Delta \rho g r^4}{\mu_p}


        Where :math:`\rho` is the density, :math:`g` is the acceleration due to gravity,
        :math:`r` is the plume conduit radius and :math:`\mu_p` is the viscosity of the plume.

        The melt flux is then obtained from:

        .. math::

           Q_m = Q_v \times F


        Where :math:`F` is the total melt fraction obtained.
        """
        Qv = _np.pi / 8 * (relative_density * 9.81 * radius**4) / viscosity
        Qm = Qv * self.F.max()
        return Qm

    def _homogenise_chemistry(self):
        """
        Homogenises the melt compositions.
        """
        # Find the elements which all melting lithologies have:
        first_lithology = True
        species = []
        for lith in self.mantle.names:
            if first_lithology is False and self.lithology_contributions[lith] > 0:
                specieslith = list(self.lithologies[lith].columns)[3:]
                species = [s for s in species if s in specieslith]
            elif first_lithology is True and self.lithology_contributions[lith] > 0:
                species = list(self.lithologies[lith].columns)[3:]
                first_lithology = False

        cm = _np.zeros([len(species)])
        for lith in self.mantle.names:
            if self.lithology_contributions[lith] > 0:
                # Get the chemistry for the lithology
                c = self.lithologies[lith][species]
                # Convert it a numpy array
                c = c.to_numpy()
                # Change nans to 0.0 to avoid affecting summation
                c = _np.nan_to_num(c, nan=0.0)
                # Get the melt fractions for the lithology
                f = self.lithologies[lith].F.to_numpy()
                # If instantaneous melts, need to pool over columns first
                for j in range(len(species)):
                    if self.MeltingColumn._species_calc_type[lith][j] == 'instantaneous':
                        cnormed = _np.zeros(_np.shape(c)[0])
                        df = f[1:] - f[: - 1]
                        for i in range(_np.shape(c)[0] - 1):
                            if f[i] > 0 and f[i - 1] > 0:
                                cnormed[i + 1] = (_np.sum(0.5 * (c[1:i + 1, j] + c[0:i, j])
                                                     * df[:i], axis=0) / f[i])
                            elif f[i] > 0:
                                cnormed[i + 1] = (_np.sum(c[1:i + 1, j]
                                                     * df[:i], axis=0) / f[i])
                            else:
                                cnormed[i + 1] = 0
                        c[:, j] = cnormed
                # Normalise melts for all lithologies
                cm += c[-1, :] * self.lithology_contributions[lith]
            self.chemistry = _pd.Series(cm, species)

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
            MeltStorageP = self.P_lithosphere
        if ShallowMeltP is None:
            ShallowMeltP = self.P_lithosphere

        T_crystallisation = self.T - (self.P - MeltStorageP) * liqdTdP

        DeepMeltTcrys = {}
        for l in self.mantle.names:
            if (self.lithologies[l].F > 0).any():
                FirstMeltRow = self.lithologies[l].F[self.lithologies[l].F > 0].idxmin()
                DeepMeltTcrys[l] = T_crystallisation.iloc[FirstMeltRow]
            else:
                DeepMeltTcrys[l] = _np.nan
        TcrysMax = _np.nanmax(list(DeepMeltTcrys.values()))

        LastMeltRow = self.P[self.P > ShallowMeltP].idxmin()
        TcrysMin = T_crystallisation.iloc[LastMeltRow]

        return TcrysMin, TcrysMax
        return TcrysMin, TcrysMax
