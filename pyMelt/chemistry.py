"""
=========
Chemistry
=========

The chemistry module provides the base classes for defining chemical elements/species for
inclusion in pyMelt calculations, alongside default implementations.
"""

import numpy as _np
import pandas as _pd
from warnings import warn as _warn

default_methods = {
    "Rb": "continuous_instantaneous",
    "Ba": "continuous_instantaneous",
    "Th": "invmel",
    "U": "invmel",
    "Nb": "invmel",
    "Ta": "invmel",
    "La": "invmel",
    "Ce": "invmel",
    "Pb": "invmel",
    "Pr": "invmel",
    "Nd": "invmel",
    "Sr": "invmel",
    "Zr": "invmel",
    "Hf": "invmel",
    "Sm": "invmel",
    "Eu": "invmel",
    "Ti": "invmel",
    "Gd": "invmel",
    "Tb": "invmel",
    "Dy": "invmel",
    "Ho": "invmel",
    "Y": "invmel",
    "Er": "invmel",
    "Yb": "invmel",
    "Lu": "invmel",
}
"""
The method to use for each element when a method isn't otherwise specified.
"""

workman05_dmm = {
    "Rb": 0.05,
    "Ba": 0.563,
    "Th": 0.0079,
    "U": 0.0032,
    "Nb": 0.1485,
    "Ta": 0.0096,
    "La": 0.192,
    "Ce": 0.550,
    "Pb": 0.018,
    "Pr": 0.107,
    "Nd": 0.581,
    "Sr": 7.664,
    "Zr": 5.082,
    "Hf": 0.157,
    "Sm": 0.239,
    "Eu": 0.096,
    "Ti": 716.3,
    "Gd": 0.358,
    "Tb": 0.070,
    "Dy": 0.505,
    "Ho": 0.115,
    "Y": 3.328,
    "Er": 0.348,
    "Yb": 0.365,
    "Lu": 0.058,
}
"""
The trace element concentrations in the depleted MORB mantle from Workman & Hart (2005). All
concentrations are in ppmw.
"""

workman05_D = {
    "Rb": 1e-5,
    "Ba": 0.00012,
    "Th": 0.001,
    "U": 0.0011,
    "Nb": 0.0034,
    "Ta": 0.0034,
    "La": 0.01,
    "Ce": 0.022,
    "Pb": 0.014,
    "Pr": 0.027,
    "Nd": 0.031,
    "Sr": 0.025,
    "Zr": 0.033,
    "Hf": 0.035,
    "Sm": 0.045,
    "Eu": 0.050,
    "Ti": 0.058,
    "Gd": 0.056,
    "Tb": 0.068,
    "Dy": 0.079,
    "Ho": 0.084,
    "Y": 0.088,
    "Er": 0.097,
    "Yb": 0.115,
    "Lu": 0.120,
}
"""
The bulk partition coefficients for MORB production from Workman & Hart (2005).
"""

stracke03_bsic = {
    "Rb": 0.57,
    "Ba": 6.59,
    "Th": 0.088,
    "U": 0.027,
    "Nb": 1.95,
    "Ta": 0.124,
    "La": 1.68,
    "Ce": 5.89,
    "Pb": 0.09,
    "Nd": 7.45,
    "Sr": 81.0,
    "Zr": 64.0,
    "Hf": 1.78,
    "Sm": 2.69,
    "Eu": 1.04,
    "Ti": 7735.0,
    "Gd": 4.03,
    "Dy": 5.01,
    "Y": 28.5,
    "Er": 3.13,
    "Yb": 2.99,
    "Lu": 0.45,
}
"""
The trace element concentrations (ppmw) in bulk subducted igneous crust from Stracke et al. (2003).
"""

palme13_pm = {
    "Rb": 0.605,
    "Ba": 6.85,
    "Th": 0.0849,
    "U": 0.0229,
    "Nb": 0.595,
    "Ta": 0.043,
    "La": 0.6832,
    "Ce": 1.7529,
    "Pb": 0.185,
    "Pr": 0.2657,
    "Nd": 1.341,
    "Sr": 22.0,
    "Zr": 10.3,
    "Hf": 0.3014,
    "Sm": 0.4347,
    "Eu": 0.1665,
    "Ti": 1265.0,
    "Gd": 0.5855,
    "Tb": 0.1075,
    "Dy": 0.7239,
    "Ho": 0.1597,
    "Y": 4.13,
    "Er": 0.4684,
    "Yb": 0.4774,
    "Lu": 0.07083,
}
"""
The composition of the primitive mantle (ppmw) from Palme & O'Neill (2013).
"""

palme13_ci = {
    "Rb": 2.32,
    "Ba": 2.42,
    "Th": 0.03,
    "U": 0.00810,
    "Nb": 0.283,
    "Ta": 0.015,
    "La": 0.2414,
    "Ce": 0.6194,
    "Pb": 2.62,
    "Pr": 0.09390,
    "Nd": 0.4737,
    "Sr": 7.79,
    "Zr": 3.63,
    "Hf": 0.1065,
    "Sm": 0.1536,
    "Eu": 0.05883,
    "Ti": 447.0,
    "Gd": 0.2069,
    "Tb": 0.03797,
    "Dy": 0.2558,
    "Ho": 0.05644,
    "Y": 1.46,
    "Er": 0.1655,
    "Yb": 0.1687,
    "Lu": 0.02503,
}
"""
Trace element concentrations in a CI chondrite (ppmw) from Palme & O'Neill (2013).
"""

# From Gibson & Geist compilation
olv_D = {
    "Rb": 0.0003,
    "Ba": 0.000005,
    "Th": 0.00005,
    "U": 0.00038,
    "Nb": 0.0005,
    "Ta": 0.0005,
    "La": 0.0005,
    "Ce": 0.0005,
    "Pb": 0.003,
    "Pr": 0.0008,
    "Nd": 0.00042,
    "Sr": 0.00004,
    "Zr": 0.0033,
    "Hf": 0.0022,
    "Sm": 0.0011,
    "Eu": 0.0016,
    "Ti": 0.015,
    "Gd": 0.0011,
    "Tb": 0.0015,
    "Dy": 0.0027,
    "Ho": 0.0016,
    "Y": 0.0099,
    "Er": 0.013,
    "Yb": 0.020,
    "Lu": 0.020,
}
"""
Trace element partition coefficients between olivine and melt, compiled by Gibson & Geist (2010).
"""

# From Gibson & Geist compilation
opx_D = {
    "Rb": 0.0002,
    "Ba": 0.000006,
    "Th": 0.002,
    "U": 0.002,
    "Nb": 0.004,
    "Ta": 0.004,
    "La": 0.0031,
    "Ce": 0.0040,
    "Pb": 0.009,
    "Pr": 0.0048,
    "Nd": 0.01200,
    "Sr": 0.0007,
    "Zr": 0.013,
    "Hf": 0.03,
    "Sm": 0.0200,
    "Eu": 0.0130,
    "Ti": 0.086,
    "Gd": 0.0130,
    "Tb": 0.0190,
    "Dy": 0.0110,
    "Ho": 0.0065,
    "Y": 0.052,
    "Er": 0.045,
    "Yb": 0.080,
    "Lu": 0.120,
}
"""
Trace element partition coefficients between orthopyroxene and melt, compiled by Gibson & Geist
(2010).
"""

# From Gibson & Geist compilation
cpx_D = {
    "Rb": 0.0004,
    "Ba": 0.0004,
    "Th": 0.0059,
    "U": 0.0094,
    "Nb": 0.015,
    "Ta": 0.015,
    "La": 0.0490,
    "Ce": 0.0800,
    "Pb": 0.012,
    "Pr": 0.126,
    "Nd": 0.17800,
    "Sr": 0.091,
    "Zr": 0.119,
    "Hf": 0.284,
    "Sm": 0.2930,
    "Eu": 0.3350,
    "Ti": 0.350,
    "Gd": 0.3500,
    "Tb": 0.4030,
    "Dy": 0.4000,
    "Ho": 0.4270,
    "Y": 0.426,
    "Er": 0.420,
    "Yb": 0.400,
    "Lu": 0.376,
}
"""
Trace element partition coefficients between clinopyroxene and melt, compiled by Gibson & Geist
(2010).
"""

# From Gibson & Geist compilation
grt_D = {
    "Rb": 0.0002,
    "Ba": 0.00007,
    "Th": 0.009,
    "U": 0.028,
    "Nb": 0.015,
    "Ta": 0.015,
    "La": 0.0010,
    "Ce": 0.0050,
    "Pb": 0.005,
    "Pr": 0.014,
    "Nd": 0.05200,
    "Sr": 0.0007,
    "Zr": 0.270,
    "Hf": 0.400,
    "Sm": 0.2500,
    "Eu": 0.4960,
    "Ti": 0.600,
    "Gd": 0.84800,
    "Tb": 1.4770,
    "Dy": 2.2000,
    "Ho": 3.3150,
    "Y": 3.100,
    "Er": 4.400,
    "Yb": 6.600,
    "Lu": 7.100,
}
"""
Trace element partition coefficients between garnet and melt, compiled by Gibson & Geist (2010).
"""

# alphaMELTS defaults
spn_D = {
    "Rb": 0.0001,
    "Ba": 0.0001,
    "Th": 0.0,
    "U": 0.0,
    "Nb": 0.0,
    "Ta": 0.0,
    "La": 0.0100,
    "Ce": 0.0100,
    "Pb": 0.0,
    "Pr": 0.01,
    "Nd": 0.0100,
    "Sr": 0.0,
    "Zr": 0.0,
    "Hf": 0.0,
    "Sm": 0.0100,
    "Eu": 0.0100,
    "Ti": 0.15,
    "Gd": 0.0100,
    "Tb": 0.0100,
    "Dy": 0.0100,
    "Ho": 0.0100,
    "Y": 0.01,
    "Er": 0.0100,
    "Yb": 0.0100,
    "Lu": 0.0100,
}
"""
Trace element partition coefficients between spinel and melt, compiled by Gibson & Geist (2010).
"""

plg_D = {
    "Rb": 0.03,
    "Ba": 0.33,
    "Th": 0.05,
    "U": 0.11,
    "Nb": 0.01,
    "Ta": 0.0,
    "La": 0.2700,
    "Ce": 0.200,
    "Pb": 0.36,
    "Pr": 0.17,
    "Nd": 0.1400,
    "Sr": 2.0,
    "Zr": 0.01,
    "Hf": 0.01,
    "Sm": 0.1100,
    "Eu": 0.7300,
    "Ti": 0.04,
    "Gd": 0.0660,
    "Tb": 0.0600,
    "Dy": 0.0550,
    "Ho": 0.0480,
    "Y": 0.03,
    "Er": 0.0100,
    "Yb": 0.031,
    "Lu": 0.0250,
}
"""
Trace element partition coefficients between plagioclase and melt, compiled by Gibson & Geist
(2010).
"""

defaultD = _pd.DataFrame(
    {
        "olv": olv_D,
        "cpx": cpx_D,
        "opx": opx_D,
        "plg": plg_D,
        "grt": grt_D,
        "spn": spn_D,
    }
)
"""
Dataframe containing the partition coefficients for each phase in INVMEL.
"""

klb1_MineralProportions = _pd.DataFrame(
    [
        [0.609, 0.125, 0.119, 0.147, 0.000, 0.000],
        [0.597, 0.233, 0.158, 0.000, 0.012, 0.000],
        [0.646, 0.208, 0.076, 0.000, 0.000, 0.070],
    ],
    columns=["olv", "opx", "cpx", "grt", "spn", "plg"],
    index=["grt_field", "spn_field", "plg_field"],
)
"""
Mineral proportions (Wt%) for KLB1 in the garnet-, spinel-, and plagioclase-field (estimated
from Jennings and Holland, 2015).
"""

kg1_MineralProportions = _pd.DataFrame(
    [
        [0.181, 0.012, 0.422, 0.385, 0.000, 0.000],
        [0.110, 0.178, 0.641, 0.000, 0.071, 0.000],
        [0.118, 0.150, 0.655, 0.000, 0.000, 0.067],
    ],
    columns=["olv", "opx", "cpx", "grt", "spn", "plg"],
    index=["grt_field", "spn_field", "plg_field"],
)
"""
Mineral proportions (Wt%) for KG1 in the garnet-, spinel-, and plagioclase-field (estimated
from Matthews et al., 2021).
"""

mo91_MineralProportions = _pd.DataFrame(
    [
        [0.598, 0.221, 0.076, 0.115, 0.000, 0.000],
        [0.578, 0.270, 0.119, 0.000, 0.033, 0.000],
        [0.636, 0.263, 0.012, 0.000, 0.000, 0.089],
    ],
    columns=["olv", "opx", "cpx", "grt", "spn", "plg"],
    index=["grt_field", "spn_field", "plg_field"],
)
"""
Mineral proportions (Wt%) for lherzolite in the garnet-, spinel-, and plagioclase-field, from
McKenzie & O'Nions (1991).
"""


class species(object):
    """
    The species base class. The class contains all the parameters and methods required to
    calculate the concentration of chemical species in melts, given T, P, F, etc.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.

    """

    def __init__(self, name, c0, **kwargs):
        self.calculation_type = None
        self.name = name
        self.c0 = c0

    def composition(self, state):
        """
        Returns the concentration of the species in the melt for the specified state.
        This function should be redefined according to the chemical model being used.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        Returns
        -------
        float
            The concentration of the species in the melt.
        """
        return _np.nan


class batchSpecies(species):
    """
    Implementation of the species class for batch melting with a constant partition coefficient.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.
    D : float
        The partition coefficient
    """

    def __init__(self, name, c0, D, **kwargs):
        self.calculation_type = "accumulated"
        self.name = name
        self.c0 = c0
        self._D = D

    def composition(self, state):
        """
        Returns the concentration in the melt during batch melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """
        return self.c0 / (self.D(state) * (1 - state["F"]) + state["F"])

    def D(self, state):
        """
        The partition coefficient. If a constant partition coefficient is used it will return that
        value. If a variable coefficient is used it will call the function to calculate it.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        Returns
        -------
        float
            The partition coefficient
        """
        if callable(self._D):
            d = self._D(state)
        else:
            d = self._D
        return d


class continuousSpecies_instantaneous(species):
    """
    Implementation of the species class for batch melting with a constant partition coefficient.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.
    D : float
        The partition coefficient
    """

    def __init__(self, name, c0, D, phi=0.005, **kwargs):
        self.calculation_type = "instantaneous"
        self.name = name
        self.c0 = c0
        self._D = D
        self.phi = phi

    def composition(self, state):
        """
        Returns the instantaneous concentration in the melt during near-fractional (continuous)
        melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """

        D = self.D(state)

        exponent = (1 - self.phi) * (1 - D) / ((1 - self.phi) * D + self.phi)

        return self.c0 / ((1 - self.phi) * D + self.phi) * (1 - state.F) ** exponent

    def D(self, state):
        """
        The partition coefficient. If a constant partition coefficient is used it will return that
        value. If a variable coefficient is used it will call the function to calculate it.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        Returns
        -------
        float
            The partition coefficient
        """
        if callable(self._D):
            d = self._D(state)
        else:
            d = self._D
        return d


class continuousSpecies_accumulated(species):
    """
    Implementation of the species class for batch melting with a constant partition coefficient.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.
    D : float
        The partition coefficient
    """

    def __init__(self, name, c0, D, phi=0.005, **kwargs):
        self.calculation_type = "accumulated"
        self.name = name
        self.c0 = c0
        self._D = D
        self.phi = phi

    def composition(self, state):
        """
        Returns the concentration in the melt during batch melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """

        Deff = (1 - self.phi) * self.D(state) + self.phi

        return self.c0 / state.F * (1 - (1 - state.F) ** (1 / Deff))

    def D(self, state):
        """
        The partition coefficient. If a constant partition coefficient is used it will return that
        value. If a variable coefficient is used it will call the function to calculate it.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        Returns
        -------
        float
            The partition coefficient
        """
        if callable(self._D):
            d = self._D(state)
        else:
            d = self._D
        return d


class mineralTransition(object):
    """
    Template for the mineralTransition class used to define mineral transitions in the invMEL
    chemistry model (and perhaps other models in the future).

    Parameters
    ----------
    parameters : dict
        A dictionary of parameters required by the transition function.
    """

    def __init__(self, parameters={}):
        self.parameters = parameters

    def transition(self, T):
        """
        Template transition method. Must be a function of temperature using parameters contained
        in the parameters dictionary.

        Parameters
        ----------
        T : float
            Temperature (in degC).

        Returns
        -------
        float
            Transition pressure (in GPa).
        """
        return _np.nan


class mineralTransition_linear(mineralTransition):
    """
    Representation of a mineral transition (e.g. garnet-in) as a linear function of temperature.

    The parameters dictionary must contain items:
        - 'intercept': the intercept pressure (T=0) in GPa
        - 'gradient': the gradient of the line (GPa/degC)
    """

    def transition(self, T):
        return self.parameters["intercept"] + self.parameters["gradient"] * T


class mineralTransition_isobaric(mineralTransition):
    """
    Representation of a mineral transition that does not depend on temperature. Temperature is
    still taken as an argument for uniformity.

    The parameters dictionary must contain items:
        - 'transition_pressure': the transition pressure in GPa
    """

    def transition(self, T):
        return self.parameters["transition_pressure"]


class invmelSpecies(species):
    """
    Implementation of the forward trace element model used by invmel (McKenzie & O'Nions, 1991).

    The default values of the mineral-melt partition coefficients are taken from REFERENCE.

    Parameters
    ----------
    name :  str
        The name of the species.
    c0 :    float
        The concentration of the species in the solid. Must be the same units as required in the
        output.
    olv_D : float
        The partition coefficient for olivine-melt.
    cpx_D : float
        The partition coefficient for clinopyroxene-melt.
    opx_D : float
        The partition coefficient for orthopyroxene-melt.
    spn_D : float
        The partition coefficient for spinel-melt.
    grt_D : float
        The partition coefficient for garnet-melt.
    plg_D : float
        The partition coefficient for plagioclase-melt.
    MineralProportions : pandas.DataFrame, default: pyMelt.chemistry.mo91_MineralProportions
        A dataframe with the proportions of each mineral phase (columns) in the garnet-, spinel-,
        and plagioclase-field for the lithology. See built in defaults for formatting of keys
        and columns: mo91_MineralProportions, klb1_MineralProportions, kg1_MineralProportions.
    density : float, default: 3.3
        The density of the mantle (g cm-3)
    cpxExhaustion : int, default: 0.18
        The melt fraction at which cpx (and grt/plg/spn) are exhausted.
    garnetOut : mineralTransition, default: 1/666.7*T + 400/666.7 GPa
        The garnet-out reaction (decompression).
    spinelIn : mineralTransition, default: 1/666.7*T + 533/666.7 GPa
        The spinel-in reaction (decompression).
    spinelOut : mineralTransition, default: 25*0.033 GPa
        The spinel-out reaction (decompression).
    plagioclaseIn : mineralTransition, default: 35*0.033 GPa
        The plagioclase-in reaction (decompression).

    Attributes
    ----------
    name : str
        Name of the species
    c0 : float
        Concentration of the species in the solid.
    D : numpy.Array
        The partition coefficients in the order: olivine, cpx, opx, spinel, garnet, plag.
    MineralProportions_solid : pandas.DataFrame
        See parameters above.
    density : float
        The density of the mantle (g cm-3)
    """

    def __init__(
        self,
        name,
        c0,
        olv_D,
        cpx_D,
        opx_D,
        spn_D,
        grt_D,
        plg_D,
        mineralProportions=mo91_MineralProportions,
        density=3.3,
        cpxExhaustion=0.18,
        garnetOut=mineralTransition_linear({'gradient': 1 / 666.7, 'intercept': 400 / 666.7}),
        spinelIn=mineralTransition_linear({'gradient': 1 / 666.7, 'intercept': 533 / 666.7}),
        spinelOut=mineralTransition_isobaric({'transition_pressure': 25.0 * 0.033}),
        plagioclaseIn=mineralTransition_isobaric({'transition_pressure': 35.0 * 0.033}),
        **kwargs
    ):
        self.calculation_type = "instantaneous"
        self.name = name
        self.c0 = c0
        self._D = {
            "olv": olv_D,
            "cpx": cpx_D,
            "opx": opx_D,
            "spn": spn_D,
            "grt": grt_D,
            "plg": plg_D,
        }
        self.mineralProportions_solid = mineralProportions
        self.density = density
        self.cpxExhaustion = cpxExhaustion
        self.garnetOut = garnetOut
        self.spinelIn = spinelIn
        self.spinelOut = spinelOut
        self.plagioclaseIn = plagioclaseIn
        self._cs = c0
        self._F_prev = 0.0
        self._cl_prev = None

    def composition(self, state):
        # Check if this is a new calculation or not:
        if state.F < self._F_prev:
            self._cs = self.c0

        if state.F == 1:
            # If the lithology is immediately fully molten:
            if self._cl_prev is None:
                return self._cs
            else:
                return self._cl_prev

        D = self.D_bulk(state["P"], state["T"], state.F)
        Pbar = self.P_bulk(state["P"], state["T"], state.F)

        if D < 1e-4:
            _warn(
                self.name
                + " is extremely incompatible, unless the step size is extremely small"
                " its partitioning behaviour is unlikely to be captured correctly. "
                "You are probably better off calculating it using the "
                "ContinuousSpecies_accumulated class."
            )

        k1 = self._dcsdX(self._F_prev, self._cs,
                         self._cl(self._cs, state.F, D, Pbar),
                         D, Pbar)
        k2 = self._dcsdX(
            self._F_prev + (state.F - self._F_prev) / 2,
            self._cs + k1 * (state.F - self._F_prev) / 2,
            self._cl(self._cs + k1 * (state.F - self._F_prev) / 2,
                     self._F_prev + (state.F - self._F_prev) / 2,
                     D, Pbar),
            D,
            Pbar,
        )
        k3 = self._dcsdX(
            self._F_prev + (state.F - self._F_prev) / 2,
            self._cs + k2 * (state.F - self._F_prev) / 2,
            self._cl(self._cs + k2 * (state.F - self._F_prev) / 2,
                     self._F_prev + (state.F - self._F_prev) / 2,
                     D, Pbar),
            D,
            Pbar,
        )
        k4 = self._dcsdX(
            self._F_prev + (state.F - self._F_prev),
            self._cs + k3 * (state.F - self._F_prev),
            self._cl(self._cs + k3 * (state.F - self._F_prev),
                     self._F_prev + (state.F - self._F_prev),
                     D, Pbar),
            D,
            Pbar,
        )
        cs = self._cs + (1 / 6) * (state.F - self._F_prev) * (k1 + 2 * k2 + 2 * k3 + k4)
        cl = self._cl(cs, state.F, D, Pbar)

        # Check if discretisation is too coarse
        if (k1 + 2 * k2 + 2 * k3 + k4) > 0 and D < 1:
            _warn(
                "Discretisation is too coarse to capture the behaviour of "
                + self.name
                + "."
            )
            cl = _np.nan
            self._cs = _np.nan

        # Prevent float errors
        elif cs < 1e-6:
            self._cs = 0
        else:
            self._cs = cs

        self._F_prev = state.F
        self._cl_prev = cl

        return cl

    def D(self, state):
        """
        The partition coefficient. If a constant partition coefficient is used it will return that
        value. If a variable coefficient is used it will call the function to calculate it.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        Returns
        -------
        float
            The partition coefficient
        """
        d = {}

        for min in self._D:
            if callable(self._D[min]):
                d[min] = self._D[min](state)
            else:
                d[min] = self._D[min]

        return d

    def mineralProportions(self, P, T):
        GarnetSpinel = self._GarnetSpinelTransition(P, T)
        SpinelPlagioclase = self._SpinelPlagioclaseTransition(P, T)
        if GarnetSpinel == 1 and SpinelPlagioclase == 1:
            mineralProportions = self.mineralProportions_solid.loc["grt_field"]
        elif GarnetSpinel == 0 and SpinelPlagioclase == 1:
            mineralProportions = self.mineralProportions_solid.loc["spn_field"]
        elif GarnetSpinel == 0 and SpinelPlagioclase == 0:
            mineralProportions = self.mineralProportions_solid.loc["plg_field"]
        elif SpinelPlagioclase == 1:
            grtField = self.mineralProportions_solid.loc["grt_field"]
            spnField = self.mineralProportions_solid.loc["spn_field"]
            mineralProportions = [
                (grtFieldProp * GarnetSpinel + spnFieldProp * (1 - GarnetSpinel))
                for grtFieldProp, spnFieldProp in zip(grtField, spnField)
            ]
            mineralProportions = _pd.Series(
                mineralProportions, index=self.mineralProportions_solid.columns
            )
        else:
            spnField = self.mineralProportions_solid.loc["spn_field"]
            plgField = self.mineralProportions_solid.loc["plg_field"]
            mineralProportions = [
                (spnFldPrp * SpinelPlagioclase + plgFldPrp * (1 - SpinelPlagioclase))
                for spnFldPrp, plgFldPrp in zip(spnField, plgField)
            ]
            mineralProportions = _pd.Series(
                mineralProportions, index=self.mineralProportions_solid.columns
            )
        return mineralProportions

    def D_bulk(self, P, T, F=0.0):
        mineralProportions = self.mineralProportions(P, T)
        Dminerals = self.D(_pd.Series({"P": P, "T": T}))

        cpxExhaustion = self.cpxExhaustion

        if F < cpxExhaustion:
            D = sum([Dminerals[min] * mineralProportions[min]
                     for min in mineralProportions.index])
        else:
            D = ((Dminerals["olv"] * mineralProportions["olv"]
                 + Dminerals["opx"] * mineralProportions["opx"])
                 / (mineralProportions["olv"] + mineralProportions["opx"]))

        return D

    def P_bulk(self, P, T, F=0.0):
        mineralProportions = self.mineralProportions(P, T)
        Dminerals = self.D(_pd.Series({"P": P, "T": T}))

        cpxExhaustion = self.cpxExhaustion

        p = {}
        if F < cpxExhaustion:
            p["cpx"] = mineralProportions["cpx"] / cpxExhaustion
            p["grt"] = mineralProportions["grt"] / cpxExhaustion
            p["spn"] = mineralProportions["spn"] / cpxExhaustion
            p["plg"] = mineralProportions["plg"] / cpxExhaustion
            p["olv"] = (
                mineralProportions["olv"]
                * (1 - (p["cpx"] + p["grt"] + p["spn"] + p["plg"]))
                / (mineralProportions["olv"] + mineralProportions["opx"])
            )
            p["opx"] = (
                mineralProportions["opx"]
                * (1 - (p["cpx"] + p["grt"] + p["spn"] + p["plg"]))
                / (mineralProportions["olv"] + mineralProportions["opx"])
            )

        else:
            p["olv"] = mineralProportions["olv"] / (
                mineralProportions["olv"] + mineralProportions["opx"]
            )
            p["opx"] = mineralProportions["opx"] / (
                mineralProportions["olv"] + mineralProportions["opx"]
            )
            p["cpx"] = 0
            p["grt"] = 0
            p["spn"] = 0
            p["plg"] = 0

        Pbar = sum([Dminerals[min] * p[min] for min in Dminerals.keys()])

        return Pbar

    def _SpinelPlagioclaseTransition(self, P, T):
        """
        Calculate proportions contributing to sp pl transition as contribution
        from the sp field mantle.

        Parameters
        ----------
        P : float
            pressure in GPa
        T : float
            temperature in degC

        Returns
        -------
        float
            contribution from sp mantle

        """
        plagIn = self.plagioclaseIn.transition(T)
        spinelOut = self.spinelOut.transition(T)

        if P <= spinelOut:
            SpinelPlagioclaseTransition = 0
        elif P >= plagIn:
            SpinelPlagioclaseTransition = 1
        else:
            SpinelPlagioclaseTransition = 1 - (plagIn - P) / (plagIn - spinelOut)
        return SpinelPlagioclaseTransition

    def _GarnetSpinelTransition(self, P, T):
        """
        Calculate proportions contributing to gt sp transition as contribution
        from the gt field mantle.

        Parameters
        ----------
        P : float
            pressure in GPa
        T : float
            temperature in degC

        Returns
        -------
        float
            contribution from gt mantle

        """
        GarnetOut = self.garnetOut.transition(T)
        SpinelIn = self.spinelIn.transition(T)
        if P >= SpinelIn:
            GarnetSpinelTransition = 1
        elif P <= GarnetOut:
            GarnetSpinelTransition = 0
        else:
            GarnetSpinelTransition = 1 - (P - SpinelIn) / (GarnetOut - SpinelIn)
        return GarnetSpinelTransition

    def _dcsdX(self, X, cs, cl, Dbar, Pbar):
        """
        Rate of change of rare-earth element concentration in point average solid.

        Parameters
        ----------
        X : float
            melt fraction
        cs : float
            point average solid residue composition
        cl : float
            liquid composition
        Dbar : float
            bulk distribution coefficient for solid assemblage
        Pbar : float
            bulk distribution coefficient for melting package

        Returns
        -------
        float
            rate of change of REE concentration in point average solid residue with respect to
            melt fraction

        """
        dcsdX = (cs - cl) / (1 - X)

        return dcsdX

    def _cl(self, cs, X, Dbar, Pbar):
        """
        Calculates instantaneous melt composition generated from a point average solid.

        Parameters
        ----------
        X : float
            melt fraction
        cs : float
            point average solid residue composition
        Dbar : float
            bulk distribution coefficient for solid assemblage
        Pbar : float
            bulk distribution coefficient for melting package

        Returns
        -------
        float
            instantaneous melt composition

        """
        cl = cs * (1 - X) / (Dbar - Pbar * X)
        return cl
