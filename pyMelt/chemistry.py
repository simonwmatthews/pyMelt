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

# From Gibson & Geist compilation
olv_D = {'Rb': 0.0003,
         'Ba': 0.000005,
         'Th': 0.00005,
         'U': 0.00038,
         'Nb': 0.0005,
         'Ta': 0.0005,
         'La': 0.0005,
         'Ce': 0.0005,
         'Pb': 0.003,
         'Pr': 0.0008,
         'Nd': 0.00042,
         'Sr': 0.00004,
         'Zr': 0.0033,
         'Hf': 0.0022,
         'Sm': 0.0011,
         'Eu': 0.0016,
         'Ti': 0.015,
         'Gd': 0.0011,
         'Tb': 0.0015,
         'Dy': 0.0027,
         'Ho': 0.0016,
         'Y': 0.0099,
         'Er': 0.013,
         'Yb': 0.020,
         'Lu': 0.020,
         }

# From Gibson & Geist compilation
opx_D = {'Rb': 0.0002,
         'Ba': 0.000006,
         'Th': 0.002,
         'U': 0.002,
         'Nb': 0.004,
         'Ta': 0.004,
         'La': 0.0031,
         'Ce': 0.0040,
         'Pb': 0.009,
         'Pr': 0.0048,
         'Nd': 0.01200,
         'Sr': 0.0007,
         'Zr': 0.013,
         'Hf': 0.03,
         'Sm': 0.0200,
         'Eu': 0.0130,
         'Ti': 0.086,
         'Gd': 0.0130,
         'Tb': 0.0190,
         'Dy': 0.0110,
         'Ho': 0.0065,
         'Y': 0.052,
         'Er': 0.045,
         'Yb': 0.080,
         'Lu': 0.120,
         }

# From Gibson & Geist compilation
cpx_D = {'Rb': 0.0004,
         'Ba': 0.0004,
         'Th': 0.0059,
         'U': 0.0094,
         'Nb': 0.015,
         'Ta': 0.015,
         'La': 0.0490,
         'Ce': 0.0800,
         'Pb': 0.012,
         'Pr': 0.126,
         'Nd': 0.17800,
         'Sr': 0.091,
         'Zr': 0.119,
         'Hf': 0.284,
         'Sm': 0.2930,
         'Eu': 0.3350,
         'Ti': 0.350,
         'Gd': 0.3500,
         'Tb': 0.4030,
         'Dy': 0.4000,
         'Ho': 0.4270,
         'Y': 0.426,
         'Er': 0.420,
         'Yb': 0.400,
         'Lu': 0.376,
         }

# From Gibson & Geist compilation
grt_D = {'Rb': 0.0002,
         'Ba': 0.00007,
         'Th': 0.009,
         'U': 0.028,
         'Nb': 0.015,
         'Ta': 0.015,
         'La': 0.0010,
         'Ce': 0.0050,
         'Pb': 0.005,
         'Pr': 0.014,
         'Nd': 0.05200,
         'Sr': 0.0007,
         'Zr': 0.270,
         'Hf': 0.400,
         'Sm': 0.2500,
         'Eu': 0.4960,
         'Ti': 0.600,
         'Gd': 0.84800,
         'Tb': 1.4770,
         'Dy': 2.2000,
         'Ho': 3.3150,
         'Y': 3.100,
         'Er': 4.400,
         'Yb': 6.600,
         'Lu': 7.100,
         }

# alphaMELTS defaults
spn_D = {'Rb': 0.0001,
         'Ba': 0.0001,
         'Th': 0.0,
         'U': 0.0,
         'Nb': 0.0,
         'Ta': 0.0,
         'La': 0.0100,
         'Ce': 0.0100,
         'Pb': 0.0,
         'Pr': 0.01,
         'Nd': 0.0100,
         'Sr': 0.0,
         'Zr': 0.0,
         'Hf': 0.0,
         'Sm': 0.0100,
         'Eu': 0.0100,
         'Ti': 0.15,
         'Gd': 0.0100,
         'Tb': 0.0100,
         'Dy': 0.0100,
         'Ho': 0.0100,
         'Y': 0.01,
         'Er': 0.0100,
         'Yb': 0.0100,
         'Lu': 0.0100,
         }

plg_D = {'Rb': 0.03,
         'Ba': 0.33,
         'Th': 0.05,
         'U': 0.11,
         'Nb': 0.01,
         'Ta': 0.0,
         'La': 0.2700,
         'Ce': 0.200,
         'Pb': 0.36,
         'Pr': 0.17,
         'Nd': 0.1400,
         'Sr': 2.0,
         'Zr':0.01,
         'Hf': 0.01,
         'Sm': 0.1100,
         'Eu': 0.7300,
         'Ti': 0.04,
         'Gd': 0.0660,
         'Tb': 0.0600,
         'Dy': 0.0550,
         'Ho': 0.0480,
         'Y': 0.03,
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

    def __init__(self, name, c0, olv_D, cpx_D, opx_D, spn_D, grt_D, plg_D,
             mineralProportions=mo91_MineralProportions, density=3.3,
             modal='NonModalVariable', modalValue=0.18,
             **kwargs):
                self.name = name
                self.c0 = c0
                self.D = {'olv':olv_D, 'cpx':cpx_D, 'opx':opx_D,
                          'spn':spn_D, 'grt':grt_D, 'plg':plg_D}
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
        GarnetSpinel = self._GarnetSpinelTransition(state.Pressure, state['T'])
        SpinelPlagioclase = self._SpinelPlagioclaseTransition(state.Pressure)
        if GarnetSpinel == 1 and SpinelPlagioclase == 1:
            mineralProportions = self.mineralProportions.loc['grt_field']
        elif GarnetSpinel == 0 and SpinelPlagioclase == 1:
            mineralProportions = self.mineralProportions.loc['spn_field']
        elif GarnetSpinel == 0 and SpinelPlagioclase == 0:
            mineralProportions = self.mineralProportions.loc['plg_field']
        elif SpinelPlagioclase == 1:
            grtField = self.mineralProportions.loc['grt_field']
            spnField = self.mineralProportions.loc['spn_field']
            mineralProportions = [(grtFieldProp * GarnetSpinel + spnFieldProp * (1 - GarnetSpinel))
                                  for grtFieldProp, spnFieldProp in zip(grtField, spnField)]
            mineralProportions = _pd.Series(mineralProportions,
                                            index=self.mineralProportions.columns)
        else:
            spnField = self.mineralProportions.loc['spn_field']
            plgField = self.mineralProportions.loc['plg_field']
            mineralProportions = [(spnFldPrp * SpinelPlagioclase
                                   + plgFldPrp * (1 - SpinelPlagioclase))
                                  for spnFldPrp, plgFldPrp in zip(spnField, plgField)]
            mineralProportions = _pd.Series(mineralProportions,
                                            index=self.mineralProportions.columns)

        modalValue = self.modalValue
        if self.modal == 'NonModalVariable':
            modalValue = (mineralProportions['cpx'] + mineralProportions['grt']
                          + mineralProportions['spn'] + mineralProportions['plg'])

        P = {}
        if state.F < modalValue:
            P['cpx'] = mineralProportions['cpx'] * state.F / modalValue
            P['grt'] = mineralProportions['grt'] * state.F / modalValue
            P['spn'] = mineralProportions['spn'] * state.F / modalValue
            P['plg'] = mineralProportions['plg'] * state.F / modalValue
            P['olv'] = (mineralProportions['olv'] * (1 - (P['cpx'] + P['grt'] + P['spn'] +
                                                     P['plg']))
                        / (mineralProportions['olv'] + mineralProportions['opx']))
            P['opx'] = (mineralProportions['opx'] * (1 - (P['cpx'] + P['grt'] + P['spn'] +
                                                     P['plg']))
                        / (mineralProportions['olv'] + mineralProportions['opx']))

            # D = sum([d * p for d, p in zip(self.D, mineralProportions)])
            D = sum([self.D[min] * mineralProportions[min] for min in mineralProportions.index])

        elif state.F >= modalValue:
            P['olv'] = (mineralProportions['olv']
                    / (mineralProportions['olv'] + mineralProportions['opx']))
            P['opx'] = (mineralProportions['opx']
                    / (mineralProportions['olv'] + mineralProportions['opx']))
            P['cpx'] = 0
            P['grt'] = 0
            P['spn'] = 0
            P['plg'] = 0
            D = ((self.D['olv'] * mineralProportions['olv']
                  + self.D['opx'] * mineralProportions['opx'])
                 / (mineralProportions['olv'] + mineralProportions['opx']))


        Pbar = sum([self.D[min] * P[min] for min in self.D.keys()])

        k1 = self._dcsdX(self._F_prev, self._cs, D, Pbar)
        k2 = self._dcsdX(self._F_prev + (state.F - self._F_prev) / 2,
                   self._cs + k1 * (state.F - self._F_prev) / 2, D, Pbar)
        k3 = self._dcsdX(self._F_prev + (state.F - self._F_prev) / 2,
                   self._cs + k2 * (state.F - self._F_prev) / 2, D, Pbar)
        k4 = self._dcsdX(self._F_prev + (state.F - self._F_prev),
                   self._cs + k3 * (state.F - self._F_prev), D, Pbar)
        self._cs = self._cs + (1 / 6) * (state.F - self._F_prev) * (k1 + 2 * k2 + 2 * k3 + k4)
        if self._cs < 1e-6:
            self._cs = 0

        cl = self._cl(self._cs, state.F, D, Pbar)

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
            SpinelPlagioclaseTransition = 1 - (35 - d) / 10
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

    def _SpinelOut(self, P):
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
        GarnetIn = self._GarnetIn(P)
        SpinelOut = self._SpinelOut(P)
        if T <= SpinelOut:
            GarnetSpinelTransition = 1
        elif T >= GarnetIn:
            GarnetSpinelTransition = 0
        else:
            GarnetSpinelTransition = 1 - (T - SpinelOut) / (GarnetIn - SpinelOut)
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
