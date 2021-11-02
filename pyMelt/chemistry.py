"""
=========
Chemistry
=========

The chemistry module provides the base classes for defining chemical elements/species for
inclusion in pyMelt calculations, alongside default implementations.
"""

import numpy as _np
import pandas as _pd

workman05_ddm = {'Rb': 0.05,
                 'Ba': 0.563,
                 'Th': 0.0079,
                 'U': 0.0032,
                 'Nb': 0.1485,
                 'Ta': 0.0096,
                 'La': 0.192,
                 'Ce': 0.550,
                 'Pb': 0.018,
                 'Pr': 0.107,
                 'Nd': 0.581,
                 'Sr': 7.664,
                 'Zr': 5.082,
                 'Hf': 0.157,
                 'Sm': 0.239,
                 'Eu': 0.096,
                 'Ti': 716.3,
                 'Gd': 0.358,
                 'Tb': 0.070,
                 'Dy': 0.505,
                 'Ho': 0.115,
                 'Y': 3.328,
                 'Er': 0.348,
                 'Yb': 0.365,
                 'Lu': 0.058
                 }

workman05_D = {'Rb': 1e-5,
               'Ba': 0.00012,
               'Th': 0.001,
               'U': 0.0011,
               'Nb': 0.0034,
               'Ta': 0.0034,
               'La': 0.01,
               'Ce': 0.022,
               'Pb': 0.014,
               'Pr': 0.027,
               'Nd': 0.031,
               'Sr': 0.025,
               'Zr': 0.033,
               'Hf': 0.035,
               'Sm': 0.045,
               'Eu': 0.050,
               'Ti': 0.058,
               'Gd': 0.056,
               'Tb': 0.068,
               'Dy': 0.079,
               'Ho': 0.084,
               'Y': 0.088,
               'Er': 0.097,
               'Yb': 0.115,
               'Lu': 0.120
               }

stracke03_bsic = {'Rb': 0.57,
                  'Ba': 6.59,
                  'Th': 0.088,
                  'U': 0.027,
                  'Nb': 1.95,
                  'Ta': 0.124,
                  'La': 1.68,
                  'Ce': 5.89,
                  'Pb': 0.09,
                  'Pr': _np.nan,
                  'Nd': 7.45,
                  'Sr': 81.0,
                  'Zr': 64.0,
                  'Hf': 1.78,
                  'Sm': 2.69,
                  'Eu': 1.04,
                  'Ti': 7735.0,
                  'Gd': 4.03,
                  'Tb': _np.nan,
                  'Dy': 5.01,
                  'Ho': _np.nan,
                  'Y': 28.5,
                  'Er': 3.13,
                  'Yb': 2.99,
                  'Lu': 0.45}

olv_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0005,
         'Ce': 0.0005,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.00042,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.0011,
         'Eu': 0.0016,
         'Ti': _np.nan,
         'Gd': 0.0011,
         'Tb': 0.0015,
         'Dy': 0.0027,
         'Ho': 0.0016,
         'Y': _np.nan,
         'Er': 0.013,
         'Yb': 0.020,
         'Lu': 0.020,
         }

opx_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0031,
         'Ce': 0.0040,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.01200,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.0200,
         'Eu': 0.0130,
         'Ti': _np.nan,
         'Gd': 0.0130,
         'Tb': 0.0190,
         'Dy': 0.0110,
         'Ho': 0.0065,
         'Y': _np.nan,
         'Er': 0.045,
         'Yb': 0.080,
         'Lu': 0.120,
         }

cpx_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0490,
         'Ce': 0.0800,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.17800,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.2930,
         'Eu': 0.3350,
         'Ti': _np.nan,
         'Gd': 0.3500,
         'Tb': 0.4030,
         'Dy': 0.4000,
         'Ho': 0.4270,
         'Y': _np.nan,
         'Er': 0.420,
         'Yb': 0.400,
         'Lu': 0.376,
         }

grt_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0010,
         'Ce': 0.0050,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.05200,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.2500,
         'Eu': 0.4960,
         'Ti': _np.nan,
         'Gd': 0.84800,
         'Tb': 1.4770,
         'Dy': 2.2000,
         'Ho': 3.3150,
         'Y': _np.nan,
         'Er': 4.400,
         'Yb': 6.600,
         'Lu': 7.100,
         }

spn_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.0100,
         'Ce': 0.0100,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.0100,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.0100,
         'Eu': 0.0100,
         'Ti': _np.nan,
         'Gd': 0.0100,
         'Tb': 0.0100,
         'Dy': 0.0100,
         'Ho': 0.0100,
         'Y': _np.nan,
         'Er': 0.0100,
         'Yb': 0.0100,
         'Lu': 0.0100,
         }

plg_D = {'Rb': _np.nan,
         'Ba': _np.nan,
         'Th': _np.nan,
         'U': _np.nan,
         'Nb': _np.nan,
         'Ta': _np.nan,
         'La': 0.2700,
         'Ce': 0.200,
         'Pb': _np.nan,
         'Pr': _np.nan,
         'Nd': 0.1400,
         'Sr': _np.nan,
         'Zr': _np.nan,
         'Hf': _np.nan,
         'Sm': 0.1100,
         'Eu': 0.7300,
         'Ti': _np.nan,
         'Gd': 0.0660,
         'Tb': 0.0600,
         'Dy': 0.0550,
         'Ho': 0.0480,
         'Y': _np.nan,
         'Er': 0.0100,
         'Yb': 0.031,
         'Lu': 0.0250,
         }

klb1_MineralProportions = _pd.DataFrame([[0.609, 0.125, 0.119, 0.147, 0.000, 0.000],
                                         [0.597, 0.233, 0.158, 0.000, 0.012, 0.000],
                                         [0.646, 0.208, 0.076, 0.000, 0.000, 0.070]],
                                        columns=['olv', 'opx', 'cpx', 'grt', 'spn', 'plg'],
                                        index=['grt_field', 'spn_field', 'plg_field'])

kg1_MineralProportions = _pd.DataFrame([[0.181, 0.012, 0.422, 0.385, 0.000, 0.000],
                                        [0.110, 0.178, 0.641, 0.000, 0.071, 0.000],
                                        [0.118, 0.150, 0.655, 0.000, 0.000, 0.067]],
                                       columns=['olv', 'opx', 'cpx', 'grt', 'spn', 'plg'],
                                       index=['grt_field', 'spn_field', 'plg_field'])

mo91_MineralProportions = _pd.DataFrame([[0.598, 0.221, 0.076, 0.115, 0.000, 0.000],
                                         [0.578, 0.270, 0.119, 0.000, 0.033, 0.000],
                                         [0.636, 0.263, 0.012, 0.000, 0.000, 0.089]],
                                        columns=['olv', 'opx', 'cpx', 'grt', 'spn', 'plg'],
                                        index=['grt_field', 'spn_field', 'plg_field'])


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
        self.name = name
        self.c0 = c0

    def instantaneous_melt(self, state):
        """
        Returns the concentration of the species in the instantaneous melt for the specified state.
        This function should be redefined according to the chemical model being used. If the model
        does not predict instantaneous melts it need not be redefined.

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

    def accumulated_melt(self, state):
        """
        Returns the concentration of the species in the accumulated melt at the specified state.
        This function should be redefined according to the chemical model being used. If the model
        does not predict accumulated melts it need not be redefined.

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


class BatchSpecies(species):
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
        self.name = name
        self.c0 = c0
        self.D = D

    def accumulated_melt(self, state):
        """
        Returns the concentration in the melt during batch melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """
        return self.c0 / (self.D * (1 - state['F']) + state['F'])


class ContinuousSpecies(species):
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
        self.name = name
        self.c0 = c0
        self.D = D
        self.phi = phi

    def instantaneous_melt(self, state):
        """
        Returns the instantaneous concentration in the melt during near-fractional (continuous)
        melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """

        exponent = (1 - self.phi) * (1 - self.D) / ((1 - self.phi) * self.D + self.phi)

        return self.c0 / ((1 - self.phi) * self.D + self.phi) * (1 - state.F)**exponent

    def accumulated_melt(self, state):
        """
        Returns the concentration in the melt during batch melting.

        Parameters
        ----------
        state : pandas.Series
            The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
            will most likely be generated automatically by the `MeltingColumn_1D` class.

        """

        Deff = (1 - self.phi) * self.D + self.phi

        return self.c0 / state.F * (1 - (1 - state.F)**(1 / Deff))


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
    densiy : float, default: 3.3
        The density of the mantle (g cm-3)
    modal : str, default: 'NonModalVariable'
        One of 'NonModalFixed', or 'NonModalVariable'. NEED MORE INFO HERE.
    modalValue : int, default: 0.18
        NEED MORE INFO HERE. The modal value to use with NonModalFixed.

    Attributes
    ----------
    name : str
        Name of the species
    c0 : float
        Concentration of the species in the solid.
    D : numpy.Array
        The partition coefficients in the order: olivine, cpx, opx, spinel, garnet, plag.
    MineralProportions : pandas.DataFrame
        See parameters above.
    density : float
        The density of the mantle (g cm-3)
    """

    def init(self, name, c0, olv_D, cpx_D, opx_D, spn_D, grt_D, plg_D,
             mineralProportions=mo91_MineralProportions, density=3.3,
             modal='NonModalVariable', modalValue=0.18,
             **kwargs):
                self.name = name
                self.c0 = c0
                self.D = np.array([olv_d, cpx_D, opx_D, spn_D, grt_D, plg_D])
                self.mineralProportions = mineralProportions
                self.density = density
                self.modal = modal
                self.modalValue = modalValue
                self._cs = c0
                self._F_prev = 0.0

    def instantaneous_melt(self, state):
        # Check if this is a new calculation or not:
        if state.F < self._F_prev:
            self._cs = c0

        if state.F == 1:
            return self._cl_prev

        # Determine the current mineralogy
        GarnetSpinel = self._GarnetSpinelTransition(state.Pressure, state.T)
        SpinelPlagioclase = self._SpinelPlagioclaseTransition(state.Pressure)
        if GarnetSpinel == 1 and SpinelPlagioclase == 1:
            mineralProportions = self.mineralProportions.loc['grt_field']
        elif GarnetSpinel == 0 and SpinelPlagioclase == 1:
            mineralProportions = self.mineralProportions.loc['spn_field']
        elif GarnetSpinel == 0 and SpinelPlagioclase == 0:
            mineralProportions = self.mineralProportions.loc['plg_field']
        elif SpinelPlagioclas == 1:
            grtField = self.mineralProportions.loc['grt_field']
            spnField = self.mineralProportions.loc['spn_field']
            mineralProportions = [(grtFieldProp * GarnetSpinel + spnFieldProp * (1 - GarnetSpinel))
                                  for grtFieldProp, spnFieldProp in zip(grtField, spnField)]
            mineralProportions = _pd.Series(mineralProportions, index=self.mineralProprtions.index)
        else:
            spnField = self.mineralProportions.loc['spn_field']
            plgField = self.mineralProportions.loc['plg_field']
            mineralProportions = [(spnFldPrp * SpinelPlagioclase
                                   + plgFldPrp * (1 - SpinelPlagioclase))
                                  for spnFldPrp, plgFldPrp in zip(SpinelField, PlagioclaseField)]
            mineralProportions = _pd.Series(mineralProportions, index=self.mineralProprtions.index)

        modalValue = self.modalValue
        if modal == 'NonModalVariable':
            modalValue = (mineralProportions['cpx'] + mineralProportions['grt']
                          + mineralProportions['spn'] + mineralProportions['plg'])

        P = np.zeros(6)
        if state.F < ModalValue:
            P[2] = mineralProportions['cpx'] * state.F / modalValue
            P[3] = mineralProportions['grt'] * state.F / modalValue
            P[4] = mineralProportions['spn'] * state.F / modalValue
            P[5] = mineralProportions['plg'] * state.F / modalValue
            P[0] = (mineralProportions['olv'] * (1 - (P[2] + P[3] + P[4] + P[5]))
                    / (mineralProportions['olv'] + mineralProportions['opx']))
            P[0] = (mineralProportions['opx'] * (1 - (P[2] + P[3] + P[4] + P[5]))
                    / (mineralProportions['olv'] + mineralProportions['opx']))

            D = sum([d * p for d, p in zip(self.D, mineralProportions)])

        elif state.F >= modalValue:
            P[0] = (mineralProportions['olv']
                    / (mineralProportions['olv'] + mineralProportions['opx']))
            P[1] = (mineralProportions['opx']
                    / (mineralProportions['olv'] + mineralProportions['opx']))
            P[2] = 0
            P[3] = 0
            P[4] = 0
            P[5] = 0
            D = ((self.D[0] * mineralProportions['olv'] + self.D[1] * mineralProportions['opx'])
                 / (mineralProportions['olv'] + mineralProportions['opx']))

        Pbar = sum([d * p for d, p in zip(self.D, P)])

        k1 = self._dcsdX(self._F_prev, self._cs, Dbar, Pbar)
        k2 = self._dcsdX(self._F_prev + (state.F - self._F_prev) / 2,
                   self._cs + k1 * (state.F - self._F_prev) / 2, Dbar, Pbar)
        k3 = self._dcsdX(self._F_prev + (state.F - self._F_prev) / 2,
                   self._cs + k2 * (state.F - self._F_prev) / 2, Dbar, Pbar)
        k4 = self._dcsdX(self._F_prev + (state.F - self._F_prev),
                   self._cs + k3 * (state.F - self._F_prev), Dbar, Pbar)
        self._cs = self._cs + (1 / 6) * (state.F - self._F_prev) * (k1 + 2 * k2 + 2 * k3 + k4)
        if self._cs < 1e-6:
            self._cs = 0.0

        cl = cl(self._cs, state.F, Dbar, Pbar)

        self._F_prev = state.F
        self._cl_prev = cl

        return cl

    def accumulated_melt(self, state):
        return _np.nan

    def _SpinelPlagioclaseTransition(self, P):
        """
        Calculate proportions contributing to sp pl transition as contribution
        from the sp field mantle.

        Parameters
        ----------
        P : float
            pressure in GPa
        rho : float
            density in g/cm3

        Returns
        -------
        float
            contribution from sp mantle

        """
        d = (P / (self.density * 10)) * 1000
        if d <= 25:
            SpinelPlagioclaseTransition = 0
        elif d >= 35:
            SpinelPlagioclaseTransition = 1
        else:
            SpinelPlagioclaseTransition = (35 - d) / 10
        return SpinelPlagioclaseTransition

    def _GarnetIn(self, P):
        """
        Calculate the temperature of garnet-in at a given pressure

        Parameters
        ----------
        P : float
            pressure in GPa

        Returns
        -------
        float
            temperature in degC

        """
        T = 666.7 * P - 400
        return T

    def SpinelOut(self, P):
        """
        Calculate the temperature of spinel-out at a given temperature

        Parameters
        ----------
        P : float
            pressure in GPa

        Returns
        -------
        float
            temperature in degC

        """
        T = 666.7 * P - 533
        return T

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
        GarnetIn = self.GarnetIn(P)
        SpinelOut = self.SpinelOut(P)
        if T <= SpinelOut:
            GarnetSpinelTransition = 1
        elif T >= GarnetIn:
            GarnetSpinelTransition = 0
        else:
            GarnetSpinelTransition = (T - SpinelOut) / (GarnetIn - SpinelOut)
        return GarnetSpinelTransition

    def _dcsdX(self, X, cs, Dbar, Pbar):
        """
        Rate of change of rare-earth element concentration in point average solid.

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
            rate of change of REE concentration in point average solid residue with respect to melt fraction

        """
        dcsdX = cs * ((1 / (1 - X)) - (1 / (Dbar - Pbar * X)))
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
        _cl : float
            instantaneous melt composition

        """
        cl = cs * (1 - X) / (Dbar - Pbar * X)
        return cl
