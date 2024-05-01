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
from pyMelt.core import InputError
from pyMelt.chemistry import matzen
from pyMelt.chemistry import data
from dataclasses import dataclass, asdict

__all__ = ['matzen', 'data']

@dataclass
class default_methods:
    """
    The method to use for each element when a method isn't otherwise specified.
    """
    Rb: str = "continuous_instantaneous"
    Ba: str = "continuous_instantaneous"
    Th: str = "invmel"
    U: str = "invmel"
    Nb: str = "invmel"
    Ta: str = "invmel"
    La: str = "invmel"
    Ce: str = "invmel"
    Pb: str = "invmel"
    Pr: str = "invmel"
    Nd: str = "invmel"
    Sr: str = "invmel"
    Zr: str = "invmel"
    Hf: str = "invmel"
    Sm: str = "invmel"
    Eu: str = "invmel"
    Ti: str = "invmel"
    Gd: str = "invmel"
    Tb: str = "invmel"
    Dy: str = "invmel"
    Ho: str = "invmel"
    Y: str = "invmel"
    Er: str = "invmel"
    Yb: str = "invmel"
    Lu: str = "invmel"




defaultD = _pd.DataFrame(
    {
        "olv": asdict(data.olv_D()),
        "cpx": asdict(data.cpx_D()),
        "opx": asdict(data.opx_D()),
        "plg": asdict(data.plg_D()),
        "grt": asdict(data.grt_D()),
        "spn": asdict(data.spn_D()),
    }
)
"""
Dataframe containing the partition coefficients for each phase in INVMEL.
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
        It can also return a dict with the compositions of the liquid plus each solid
        phase included.

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
    
    def mineralProportions(self, state):
        """
        Returns the proportions of mineral phases in the solid residue for the specified
        state. This function may be redefined according to the chemical model being used.
        The simplest models (e.g., batchSpecies or fractionalSpecies) do not imply any 
        particular residue composition, so this method need not be redefined, and NaN can
        be returned.

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
        self.calculation_type = "liquidConcentrationAccumulated"
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
        self.calculation_type = "liquidConcentrationInstantaneous"
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
        self.calculation_type = "liquidConcentrationAccumulated"
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

    # The defaults are set to None so that the class can be used for calculating minineral
    # proportions without having to provide fake information.
    def __init__(
        self,
        name=None,
        c0=None,
        olv_D=None,
        cpx_D=None,
        opx_D=None,
        spn_D=None,
        grt_D=None,
        plg_D=None,
        mineralProportions=data.mo91_MineralProportions,
        density=3.3,
        cpxExhaustion=0.18,
        garnetOut=mineralTransition_linear({'gradient': 1 / 666.7, 'intercept': 400 / 666.7}),
        spinelIn=mineralTransition_linear({'gradient': 1 / 666.7, 'intercept': 533 / 666.7}),
        spinelOut=mineralTransition_isobaric({'transition_pressure': 25.0 * 0.033}),
        plagioclaseIn=mineralTransition_isobaric({'transition_pressure': 35.0 * 0.033}),
        **kwargs
    ):
        self.calculation_type = "liquidConcentrationInstantaneous"
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

        # Check if discretisation is too course
        if (k1 + 2 * k2 + 2 * k3 + k4) > 0 and D < 1:
            _warn(
                "Discretisation is too course to capture the behaviour of "
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

    def mineralProportions(self, state):
        P = state['P']
        T = state['T']
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
        mineralProportions = self.mineralProportions(_pd.Series({"P": P, "T": T}))
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
        mineralProportions = self.mineralProportions(_pd.Series({"P": P, "T": T}))
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

class phaseDiagramTraceSpecies(species):
    """
    Calculates trace element partitioning based on a phase diagram. Needs more info.
    """
    # The defaults are set to None so that the class can be used for calculating minineral
    # proportions without having to provide fake information.
    def __init__(self, name=None, c0=None, phaseDiagram=None,
                 olv_D=None, cpx_D=None, opx_D=None, spn_D=None, grt_D=None, plg_D=None,
                 porosity=0.0, **kwargs):
        self.calculation_type = "liquidConcentrationInstantaneous"
        self.name = name
        self.c0 = c0
        self.phaseDiagram = phaseDiagram
        self._D = {"olv": olv_D,
                   "cpx": cpx_D,
                   "opx": opx_D,
                   "spn": spn_D,
                   "grt": grt_D,
                   "plg": plg_D}
        self.porosity = porosity
        self._cs = c0
        self._F_prev = 0.0
        self._cl_prev = None
        self._kwargs = kwargs
        self._kwargs['phaseDiagram'] = self.phaseDiagram


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

        D = self.D_bulk(state)

        k1 = self._dcsdX(self._F_prev, self._cs,
                         self._cl(self._cs, D))
        k2 = self._dcsdX(
            self._F_prev + (state.F - self._F_prev) / 2,
            self._cs + k1 * (state.F - self._F_prev) / 2,
            self._cl(self._cs + k1 * (state.F - self._F_prev) / 2,
                     D),
        )
        k3 = self._dcsdX(
            self._F_prev + (state.F - self._F_prev) / 2,
            self._cs + k2 * (state.F - self._F_prev) / 2,
            self._cl(self._cs + k2 * (state.F - self._F_prev) / 2,
                     D),
        )
        k4 = self._dcsdX(
            self._F_prev + (state.F - self._F_prev),
            self._cs + k3 * (state.F - self._F_prev),
            self._cl(self._cs + k3 * (state.F - self._F_prev),
                     D),
        )
        cs = self._cs + (1 / 6) * (state.F - self._F_prev) * (k1 + 2 * k2 + 2 * k3 + k4)
        cl = self._cl(cs, D)


        # Check if discretisation is too course
        if (k1 + 2 * k2 + 2 * k3 + k4) > 0 and D < 1:
            _warn(
                "Discretisation is too course to capture the behaviour of "
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
        
        # Calculate mineral compositions
        c = {'liq': cl}
        for min in D:
            c[min] = cl * D[min]

        self._F_prev = state.F
        self._cl_prev = cl

        return c

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
                d[min] = self._D[min](state, mineral=min, element=self.name, **self._kwargs)
            else:
                d[min] = self._D[min]

        return d

    def D_bulk(self, state):
        """
        """
        Dminerals = self.D(state)
        mineralProportions = self.mineralProportions(state)

        D = 0.0
        for mineral in Dminerals.keys():
            if mineralProportions[mineral] > 1e-10:
                D += Dminerals[mineral] * mineralProportions[mineral]

        D = (1 - self.porosity) * D + self.porosity

        return D

    def mineralProportions(self, state):
        props = {}
        total = 0.0
        for mineral in self._D.keys():
            props[mineral] = self.phaseDiagram(mineral + '_mass', state)
            total += props[mineral] 
        for mineral in self._D.keys():
            props[mineral] = props[mineral] / total
        return props

    def _dcsdX(self, X, cs, cl):
        """
        Rate of change of rare-earth element concentration in point average solid. This
        expression is only strictly valid in the limit dX->0.

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

    def _cl(self, cs, Dbar):
        """
        Calculates instantaneous melt composition generated from a point average solid.

        Parameters
        ----------
        cs : float
            point average solid residue composition
        Dbar : float
            bulk distribution coefficient for solid assemblage


        Returns
        -------
        float
            instantaneous melt composition

        """
        cl = cs / Dbar
        return cl

class phaseDiagramMajorSpecies(species):
    """
    Calculates trace element partitioning based on a phase diagram. Needs more info.
    """
    # The defaults are set to None so that the class can be used for calculating minineral
    # proportions without having to provide fake information.
    def __init__(self, name=None, phaseDiagram=None, suffix='_wtpt',
                 **kwargs):
        self.calculation_type = "liquidConcentrationInstantaneous"
        self.name = name
        self.phaseDiagram = phaseDiagram
        self.suffix = suffix
        self._F_prev = 0.0

    def composition(self, state):
        c = {'liq': self.phaseDiagram((self.phaseDiagram.liquid_name + '_' 
                                      + self.name + self.suffix), state)}
        for mineral in self.phaseDiagram.minerals:
            try:
                c[mineral] = self.phaseDiagram(mineral + '_' + self.name + self.suffix, 
                                               state)
            except Exception:
                _warn('The ' + self.name + ' of ' + mineral + ' could not be found.')
        return c

    def mineralProportions(self, state):
        props = {}
        total = 0.0
        for mineral in self.phaseDiagram.minerals:
            props[mineral] = self.phaseDiagram(mineral + '_mass', state)
            total += props[mineral] 
        for mineral in self.phaseDiagram.minerals:
            props[mineral] = props[mineral] / total
        return props


# FUNCTIONS FOR LATTICE STRAIN CALCULATIONS

def lattice_D(state, D0, r0, ri, Em, **kwargs):
    """
    Calculates the partition coefficient based on lattice strain parameters.

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    D0 : float
        The reference partition coefficient
    r0 : float
        The reference site ionic radius (in m)
    ri : float
        The ionic radius of the element (in m)
    Em : float
        The Young's modulus of the lattice site (in Pa)

    Returns
    -------
    float
        The partition coefficient
    """
    Na = 6.02214076e23
    R = 8.31446261815324
    T = state['T'] + 273.15

    if callable(D0):
        D0 = D0(state, **kwargs)

    if callable(r0):
        r0 = r0(state, **kwargs)

    if callable(ri):
        ri = ri(state, **kwargs)

    if callable(Em):
        Em = Em(state, **kwargs)

    try:
        Di = (D0 * _np.exp(-4 * _np.pi * Na * Em *
                           (0.5 * r0 * (ri - r0)**2 + 1/3 * (ri-r0)**3)
                           / ( R * T)))
    except Exception:
        Di = _np.nan
    return Di

def lattice_cpx_D(state, cpx_D0, cpx_r0, ri, cpx_Em, **kwargs):
    """
    Redefines latticeD with parameters named specifically for cpx.

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    cpx_D0 : float
        The reference partition coefficient
    cpx_r0 : float
        The reference site ionic radius (in m)
    ri : float
        The ionic radius of the element (in m)
    cpx_Em : float
        The Young's modulus of the lattice site (in Pa)

    Returns
    -------
    float
        The partition coefficient
    """
    return lattice_D(state, cpx_D0, cpx_r0, ri, cpx_Em, **kwargs)

def lattice_grt_D(state, grt_D0, grt_r0, ri, grt_Em, **kwargs):
    """
    Redefines latticeD with parameters named specifically for garnet.

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    grt_D0 : float
        The reference partition coefficient
    grt_r0 : float
        The reference site ionic radius (in m)
    ri : float
        The ionic radius of the element (in m)
    grt_Em : float
        The Young's modulus of the lattice site (in Pa)

    Returns
    -------
    float
        The partition coefficient
    """
    return lattice_D(state, grt_D0, grt_r0, ri, grt_Em, **kwargs)

def lattice_plg_D(state, plg_D0, plg_r0, ri, plg_Em, **kwargs):
    """
    Redefines latticeD with parameters named specifically for plagioclase.

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    plg_D0 : float
        The reference partition coefficient
    plg_r0 : float
        The reference site ionic radius (in m)
    ri : float
        The ionic radius of the element (in m)
    plg_Em : float
        The Young's modulus of the lattice site (in Pa)

    Returns
    -------
    float
        The partition coefficient
    """
    return lattice_D(state, plg_D0, plg_r0, ri, plg_Em, **kwargs)

def lattice_cpx_r0(state, z, cpx_M2_Ca=None, cpx_M1_Al=None, **kwargs):
    """
    Calculates cpx reference ionic radius for +3 ions from its composition. From Wood & Blundy
    (1997).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion. Only 3 is supported currently.
    cpx_M2_Ca : float, None, or function.
        The mole fraction of Ca on the M2 site. If None, it will look for a phaseDiagram object in
        kwargs with a parameter called 'cpx_xCaM2'.
    cpx_M1_Al : float or None, or function.
        The mole fraction of Al on the M1 site. If None, it will look for a phaseDiagram object in
        kwargs with a parameter called 'cpx_xAlM1'.

    Returns
    -------
    float
        The cpx reference ionic radius (in m).
    """
    if callable(cpx_M2_Ca):
        cpx_M2_Ca = cpx_M2_Ca(state, **kwargs)
    elif cpx_M2_Ca is None and 'phaseDiagram' in kwargs:
        cpx_M2_Ca = kwargs['phaseDiagram']('cpx_xCaM2',  state)

    if callable(cpx_M1_Al):
        cpx_M1_Al = cpx_M1_Al(state, **kwargs)
    elif cpx_M1_Al is None and 'phaseDiagram' in kwargs:
        cpx_M1_Al = kwargs['phaseDiagram']('cpx_xAlM1',  state)

    if z == 3:
        r0 = (0.974 + 0.067*cpx_M2_Ca - 0.051*cpx_M1_Al) * 1e-10
    else:
        r0 = _np.nan
    return r0

def lattice_cpx_Em(state, z, **kwargs):
    """
    Calculates cpx Young's modulus for +3 ions at specified pressure and temperature.

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion. Only 3 is supported currently.

    Returns
    -------
    float
        The Young's modulus (in Pa).
    """
    if z == 3:
        Em = (318.6 + 6.9*state.P - 0.036*(state['T'] + 273.15))*1e9
    else:
        Em = _np.nan
    return Em

def lattice_cpx_D0(state, z, melt_Mgn=None, cpx_M1_Mg=None, **kwargs):
    """
    Calculates cpx reference partition coefficient for +3 ions from its composition. From Wood &
    Blundy (1997).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion. Only 3 is supported currently.
    melt_Mgn : float, None, or function.
        The Mg# of the melt in equilibrium with the crystal. If None, it will look for a
        phaseDiagram object in kwargs with a parameter called 'liq_Mg#'.
    cpx_M1_Mg : float, None, or function.
        The mole fraction of Mg on the M1 site. If None, it will look for a phaseDiagram object in
        kwargs with a parameter called 'cpx_xMgM1'.

    Returns
    -------
    float
        The reference partition coefficient.
    """
    R = 8.31446261815324

    if callable(melt_Mgn):
        melt_Mgn = melt_Mgn(state, **kwargs)
    elif melt_Mgn is None and 'phaseDiagram' in kwargs:
        melt_Mgn = kwargs['phaseDiagram']('liq_Mg#',  state)

    if callable(cpx_M1_Mg):
        cpx_M1_Mg = cpx_M1_Mg(state, **kwargs)
    elif cpx_M1_Mg is None and 'phaseDiagram' in kwargs:
        cpx_M1_Mg = kwargs['phaseDiagram']('cpx_xMgM1',  state)

    if z == 3:
        if cpx_M1_Mg > 0:
            D = (melt_Mgn / cpx_M1_Mg
                 * _np.exp((88750 - 65.644 * (state['T'] + 273.15) + 7050 * state.P
                            - 770 * (state.P)**2)
                          / (R * (state['T'] + 273.15))))
        else:
            D = 0
    else:
        D = _np.nan

    return D


def lattice_plg_D0(state, z, plg_an=None, plg_ab=None, **kwargs):
    """
    Calculates the reference partition coefficient for plagioclase according to Sun et al. (2017).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion, must be between 1 and 3.
    plg_an : float or None, default: None
        The mole fraction of anorthite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'pl_anorthite'.
    plg_ab : float or None, default: None
        The mole fraction of albite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'pl_albite'.

    Returns
    -------
    float
        The reference partition coefficient.
    """
    R = 8.31446261815324

    if callable(plg_an):
        plg_an = plg_an(state, **kwargs)
    elif plg_an is None and 'phaseDiagram' in kwargs:
        plg_an = kwargs['phaseDiagram']('pl_anorthite',  state)

    if callable(plg_ab):
        plg_ab = plg_ab(state, **kwargs)
    elif plg_ab is None and 'phaseDiagram' in kwargs:
        plg_ab = kwargs['phaseDiagram']('pl_albite',  state)

    if z == 3:
        lnD = (16.05 - ((19.45 + 1.17 * state.P**2)/(R * (state['T'] + 273.15)) * 1e4)
               - 5.17 * plg_an**2)
    elif z == 2:
        lnD = ((6910 - 2542 * state.P**2) / (R * (state['T'] + 273.15))
               + 2.39 * plg_ab**2)
    elif z == 1:
        lnD = (-9.99 + (11.37 + 0.49 * state['P']) / (R * (state['T'] + 273.15)) * 1e4
               + 1.7 * plg_an**2)
    else:
        return _np.nan

    return _np.exp(lnD)

def lattice_plg_r0(state, z, plg_ab=None, **kwargs):
    """
    Calculates the reference ionic radius (in m) in plagioclase for the specified ion charge and
    plagioclase composition, according to Sun et al. (2017).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion, must be between 1 and 3.
    plg_ab : float or None, default: None
        The mole fraction of albite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'pl_albite'.

    Returns
    -------
    float
        The reference ionic radius (in m).
    """

    if z == 3:
        r0 = 1.179
    elif z == 2:
        if callable(plg_ab):
            plg_ab = plg_ab(state, **kwargs)
        elif plg_ab is None and 'phaseDiagram' in kwargs:
            plg_ab = kwargs['phaseDiagram']('pl_albite',  state)

        r0 = 1.189 + 0.075 * plg_ab
    elif z == 1:
        r0 = 1.213
    else:
        r0 = _np.nan

    return r0 * 1e-10

def lattice_plg_Em(state, z, plg_ab=None, **kwargs):
    """
    Returns the Young's modulus (in Pa) for the plagioclase lattice sites, depending on ion charge
    and plagioclase composition, according to Sun et al. (2017).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion, must be between 1 and 3.
    plg_ab : float or None, default: None
        The mole fraction of albite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'pl_albite'.

    Returns
    -------
    float
        The Young's modulus (in Pa).
    """


    if z == 3:
        E = 196
    elif z == 2:
        if callable(plg_ab):
            plg_ab = plg_ab(state, **kwargs)
        elif plg_ab is None and 'phaseDiagram' in kwargs:
            plg_ab = kwargs['phaseDiagram']('pl_albite',  state)

        E = 719 - 487 * lattice_plg_r0(state, z, plg_ab) * 1e10
    elif z == 1:
        E = 47
    else:
        E = _np.nan

    return E * 1e9

def lattice_grt_Em(state, z, grt_prp=None, grt_grs=None, grt_alm=None, grt_sps=None, grt_and=None,
                   grt_uvr=None, **kwargs):
    """
    Calculates the Young's modulus (in Pa) for garnet, according to Westrenen & Draperc(2007).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion. Only 3 is currently supported.
    grt_prp : float, None, or function, default:  None
        The mole fraction of pyrope. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_pyrope'.
    grt_grs : float, None, or function, default:  None
        The mole fraction of grossular. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_grossular'.
    grt_alm : float, None, or function, default:  None
        The mole fraction of almandine. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_almandine'.
    grt_sps : float, None, or function, default:  None
        The mole fraction of spessartine. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_spessartine'.
    grt_and : float, None, or function, default:  None
        The mole fraction of andradite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_andradite'.
    grt_uvr : float, None, or function, default:  None
        The mole fraction of uvarovite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_uvarovite'.

    Returns
    -------
    float
        Young's modulus (in Pa).
    """
    if callable(grt_prp):
        grt_prp = grt_prp(state, **kwargs)
    elif grt_prp is None and 'phaseDiagram' in kwargs:
        grt_prp = kwargs['phaseDiagram']('g_pyrope',  state)

    if callable(grt_grs):
        grt_grs = grt_grs(state, **kwargs)
    elif grt_grs is None and 'phaseDiagram' in kwargs:
        grt_grs = kwargs['phaseDiagram']('g_grossular',  state)

    if callable(grt_alm):
        grt_alm = grt_alm(state, **kwargs)
    elif grt_alm is None and 'phaseDiagram' in kwargs:
        grt_alm = kwargs['phaseDiagram']('g_almandine',  state)

    if callable(grt_sps):
        grt_sps = grt_sps(state, **kwargs)
    elif grt_sps is None and 'phaseDiagram' in kwargs:
        grt_sps = kwargs['phaseDiagram']('g_spessartine',  state)

    if callable(grt_and):
        grt_and = grt_and(state, **kwargs)
    elif grt_and is None and 'phaseDiagram' in kwargs:
        grt_and = kwargs['phaseDiagram']('g_andradite',  state)

    if callable(grt_uvr):
        grt_uvr = grt_uvr(state, **kwargs)
    elif grt_uvr is None and 'phaseDiagram' in kwargs:
        grt_uvr = kwargs['phaseDiagram']('g_uvarovite',  state)

    if z == 3:
        r0 = lattice_grt_r0(state, z, grt_prp, grt_grs, grt_alm, grt_sps, grt_and, grt_uvr, **kwargs)
        al_pfu = 2 * (grt_prp + grt_alm + grt_sps + grt_grs)
        cr_pfu = 2 * grt_uvr
        Em = (2826 * (1.38 + r0*1e10)**-3 + 12.4 * state['P'] - 0.072 * (state['T'] + 273.15)
              + 237 * (al_pfu + cr_pfu))
    else:
        Em = _np.nan

    return Em * 1e9



def lattice_grt_r0(state, z, grt_prp=None, grt_grs=None, grt_alm=None, grt_sps=None, grt_and=None,
                   grt_uvr=None, **kwargs):
    """
    Calculates the reference ionic radius (in m) for garnet, according to Westrenen & Draper
    (2007).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion. Only 3 is currently supported.
    grt_prp : float, None, or function, default:  None
        The mole fraction of pyrope. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_pyrope'.
    grt_grs : float, None, or function, default:  None
        The mole fraction of grossular. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_grossular'.
    grt_alm : float, None, or function, default:  None
        The mole fraction of almandine. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_almandine'.
    grt_sps : float, None, or function, default:  None
        The mole fraction of spessartine. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_spessartine'.
    grt_and : float, None, or function, default:  None
        The mole fraction of andradite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_andradite'.
    grt_uvr : float, None, or function, default:  None
        The mole fraction of uvarovite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_uvarovite'.

    Returns
    -------
    float
        The reference ionic radius (in m).
    """
    if callable(grt_prp):
        grt_prp = grt_prp(state, **kwargs)
    elif grt_prp is None and 'phaseDiagram' in kwargs:
        grt_prp = kwargs['phaseDiagram']('g_pyrope',  state)

    if callable(grt_grs):
        grt_grs = grt_grs(state, **kwargs)
    elif grt_grs is None and 'phaseDiagram' in kwargs:
        grt_grs = kwargs['phaseDiagram']('g_grossular',  state)

    if callable(grt_alm):
        grt_alm = grt_alm(state, **kwargs)
    elif grt_alm is None and 'phaseDiagram' in kwargs:
        grt_alm = kwargs['phaseDiagram']('g_almandine',  state)

    if callable(grt_sps):
        grt_sps = grt_sps(state, **kwargs)
    elif grt_sps is None and 'phaseDiagram' in kwargs:
        grt_sps = kwargs['phaseDiagram']('g_spessartine',  state)

    if callable(grt_and):
        grt_and = grt_and(state, **kwargs)
    elif grt_and is None and 'phaseDiagram' in kwargs:
        grt_and = kwargs['phaseDiagram']('g_andradite',  state)

    if callable(grt_uvr):
        grt_uvr = grt_uvr(state, **kwargs)
    elif grt_uvr is None and 'phaseDiagram' in kwargs:
        grt_uvr = kwargs['phaseDiagram']('g_uvarovite',  state)

    if z == 3:
        r0 = (0.9302 * grt_prp + 0.993 * grt_grs + 0.916 * grt_alm + 0.946 * grt_sps
              + 1.05 * (grt_and + grt_uvr) - 0.0044 * (state.P - 3)
              + 0.000058 * (state['T'] + 273.15 - 1818)
              )
    else:
        r0 = _np.nan

    return r0*1e-10

def lattice_grt_D0(state, z, grt_FeO=None, liq_FeO=None, grt_grs=None, grt_and=None,
                  grt_uvr=None, **kwargs):
    """
    Calculates the reference partition coefficient for garnet, according to Westrenen & Draper
    (2007).

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    z : int
        The charge on the ion. Only 3 is currently supported.
    grt_FeO : float, None, or function, default:  None
        The wt% of FeO in garnet. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_FeO_wtpt'.
    liq_FeO : float, None, or function, default:  None
        The wt% of FeO in the melt. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'liq_FeO_wtpt'.
    grt_grs : float, None, or function, default:  None
        The mole fraction of grossular. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_grossular'.
    grt_and : float, None, or function, default:  None
        The mole fraction of andradite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_andradite'.
    grt_uvr : float, None, or function, default:  None
        The mole fraction of uvarovite. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'g_uvarovite'.

    Returns
    -------
    float
        The reference ionic radius (in m).
    """
    R = 8.31446261815324

    if callable(grt_FeO):
        grt_FeO = grt_FeO(state, **kwargs)
    elif grt_FeO is None and 'phaseDiagram' in kwargs:
        grt_FeO = kwargs['phaseDiagram']('g_FeO_wtpt',  state)

    if callable(liq_FeO):
        liq_FeO = liq_FeO(state, **kwargs)
    elif liq_FeO is None and 'phaseDiagram' in kwargs:
        liq_FeO = kwargs['phaseDiagram']('liq_FeO_wtpt',  state)

    if callable(grt_grs):
        grt_grs = grt_grs(state, **kwargs)
    elif grt_grs is None and 'phaseDiagram' in kwargs:
        grt_grs = kwargs['phaseDiagram']('g_grossular',  state)

    if callable(grt_and):
        grt_and = grt_and(state, **kwargs)
    elif grt_and is None and 'phaseDiagram' in kwargs:
        grt_and = kwargs['phaseDiagram']('g_andradite',  state)

    if callable(grt_uvr):
        grt_uvr = grt_uvr(state, **kwargs)
    elif grt_uvr is None and 'phaseDiagram' in kwargs:
        grt_uvr = kwargs['phaseDiagram']('g_uvarovite',  state)

    if z == 3:
        if grt_FeO > 0 and liq_FeO > 0:
            gamma_Fe = _np.exp((19000 * (grt_grs + grt_and + grt_uvr)**2)
                               / (R * (state['T'] + 273.15)))
            DFe = grt_FeO / liq_FeO

            D0 = (_np.exp((400290 + 4586 * state['P'] - 218 * (state['T'] + 273.15))
                          / (R * (state['T'] + 273.15)))
                  / (gamma_Fe * DFe)**2)
        else:
            D0 = 0.0
    else:
        D0 = _np.nan

    return D0

