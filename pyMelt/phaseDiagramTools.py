"""
===================
Phase Diagram Tools
===================

Contains classes and methods for importing and using phase diagrams to do chemistry calculations.

"""

import numpy as _np
from scipy.interpolate import interp2d as _interp2d
from scipy.interpolate import interp1d as _interp1d
from copy import copy as _copy
import pandas as _pd
import os
import pickle
from json import load
import matplotlib.pyplot as _plt

class gridsMelts(object):
    """
    Methods for processing MELTS input.
    """
    def __init__(self,
                 df,
                 phases=['ol', 'cpx', 'opx', 'g', 'spn', 'liq', 'pl'],
                 oxides={'SiO2':  28.085 + 15.999 * 2,
                         'MgO':   24.305 + 15.999,
                         'FeO':   55.845 + 15.999,
                         'CaO':   40.078 + 15.999,
                         'Al2O3': 2 * 26.982 + 15.999 * 3,
                         'Na2O':  22.99 * 2 + 15.999,
                         'K2O':   39.098 * 2 + 15.999,
                         'MnO':   54.938 + 15.999,
                         'TiO2':  79.867,
                         'P2O5':  2 * 30.974 + 5 * 15.999,
                         'Cr2O3': 151.992,
                         'NiO':   58.693 + 16,
                         'CoO':   44.01,
                         'Fe2O3': 55.845 * 2 + 15.999 * 3,
                         'H2O':   18.02,
                         'CO2':   44.01,
                         'O': 15.999}):

        self.df = df.copy()
        self.phases  = phases
        self.oxides = oxides


    def make_f_grid(self, variable, deltaf=0.01, fill_value=[0.0, 0.0], use_edge_values=False):
        f = _np.arange(0.0, 1.0, deltaf)
        p = self.df['pressure'].unique()
        pp, ff = _np.meshgrid(p, f)
        grid = _np.zeros([3, _np.shape(pp)[0], _np.shape(pp)[1]])
        grid[0, :, :] = pp
        grid[1, :, :] = ff
        fill_value = _copy(fill_value)

        for i in range(_np.shape(pp)[1]):
            dfp = self.df[(self.df.pressure==pp[0, i])]
            x = dfp.liq_mass.copy()
            y = dfp[variable].copy()

            # Check that there's a liq=0 value, if not create one:
            if _np.isnan(x.iloc[0]) == False and x.iloc[0] > 0:
                ind = list(x.index)
                ind = [ind[0]-1] + ind
                xnew = _np.zeros(_np.shape(x)[0]+1)
                xnew[1:] = x
                xnew[0] = _np.nan
                x = _pd.Series(xnew, index=ind)
                ynew = _np.zeros(_np.shape(y)[0]+1)
                ynew[1:] = y
                ynew[0] = _np.nan
                y = _pd.Series(ynew, index=ind)


            # To ensure smoothness to zero, we will artificially overshoot F=0:
            first_finite_ind = x[(_np.isnan(x) == False)].index[0]
            dfdi0 = x.loc[first_finite_ind + 1] - x.loc[first_finite_ind]
            # To accomodate phase diagrams where the first increments of melting are not present:
            dfdi_multiplier = 1
            while x.loc[first_finite_ind] - dfdi0 * dfdi_multiplier > 0:
                dfdi_multiplier += 1
            x.loc[first_finite_ind - 1] = x.loc[first_finite_ind] - dfdi0 * dfdi_multiplier

            # Now ensure that the variable can overshoot on both ends:
            y = y[(_np.isnan(x) == False)]

            if len(y[(_np.isnan(y) == False)]) > 1:
                # Check lower end:
                if _np.isnan(y.iloc[0]):
                    first_finite_ind = y[(_np.isnan(y) == False)].index[0]
                    dydx = (y.loc[first_finite_ind + 1] - y.loc[first_finite_ind])/(x.loc[first_finite_ind + 1] - x.loc[first_finite_ind])
                    y.loc[first_finite_ind - 1] = y.loc[first_finite_ind] + dydx*(x.loc[first_finite_ind - 1] - x.loc[first_finite_ind]) * dfdi_multiplier
                    if use_edge_values is True:
                        fill_value[0] = y.loc[first_finite_ind - 1]

                # Check upper end:
                if _np.isnan(y.iloc[len(y)-1]):
                    last_finite_ind = y[(_np.isnan(y) == False)].index[-1]
                    dydx = (y.loc[last_finite_ind] - y.loc[last_finite_ind - 1])/(x.loc[last_finite_ind] - x.loc[last_finite_ind - 1])
                    y.loc[last_finite_ind + 1] = y.loc[last_finite_ind] + dydx*(x.loc[last_finite_ind + 1] - x.loc[last_finite_ind]) * dfdi_multiplier
                    if use_edge_values is True:
                        fill_value[1] = y.loc[last_finite_ind + 1]

                # The fill value must be a tuple for interp1d:
                if isinstance(fill_value, tuple):
                    fill_value_use = fill_value
                else:
                    fill_value_use = tuple(fill_value)

                x = x[(_np.isnan(x) == False) & (_np.isnan(y) == False)]

                if len(x) > 1:
                    y = y[(_np.isnan(x) == False) & (_np.isnan(y) == False)]

                    try:
                        pinterp = _interp1d(x, y, fill_value=fill_value_use, bounds_error=False, kind='quintic')
                    except:
                        try:
                            pinterp = _interp1d(x, y, fill_value=fill_value_use, bounds_error=False, kind='cubic')
                        except:
                            pinterp = _interp1d(x, y, fill_value=fill_value_use, bounds_error=False)

                    grid[2, :, i] = pinterp(f)

        return grid

class gridsThermocalc(object):
    """
    Methods for processing thermocalc input.
    """
    def __init__(self,
                 df,
                 phases=['olv', 'cpx', 'opx', 'grt', 'spn', 'liq', 'plg'],
                 oxides={'SiO2':  28.085 + 15.999 * 2,
                         'MgO':   24.305 + 15.999,
                         'FeO':   55.845 + 15.999,
                         'CaO':   40.078 + 15.999,
                         'Al2O3': 2 * 26.982 + 15.999 * 3,
                         'Na2O':  22.99 * 2 + 15.999,
                         'K2O':   39.098 * 2 + 15.999,
                         'MnO':   54.938 + 15.999,
                         'TiO2':  79.867,
                         'P2O5':  2 * 30.974 + 5 * 15.999,
                         'Cr2O3': 151.992,
                         'NiO':   58.693 + 16,
                         'CoO':   44.01,
                         'Fe2O3': 55.845 * 2 + 15.999 * 3,
                         'H2O':   18.02,
                         'CO2':   44.01,
                         'O': 15.999}):

        self.df = df.copy()
        self.phases  = phases
        self.oxides = oxides

        # Copying the dataframe will prevent fragmentation warnings and improve performance.
        self.calculate_mass_fractions()
        self.df = self.df.copy()
        self.calculate_phase_wtpt()
        self.df = self.df.copy()
        self.calculate_phase_Mgn()
        self.df = self.df.copy()

        # Thermocalc uses kbar as the pressure unit. pyMelt uses GPa, convert:
        self.df['pressure'] = self.df['pressure'] / 10

    def calculate_mass_fractions(self):
        """
        """

        # TC calculates modes on a 1-atom basis.

        # First, calculate un-normalised mass fractions
        phases_found = []
        for ph in self.phases:
            if ph in self.df.columns:
                phases_found.append(ph)
                mass = _np.zeros([_np.shape(self.df)[0]])
                for ox in self.oxides.keys():
                    if ph + '_' + ox in self.df.columns:
                        mass += self.df[ph + '_' + ox] * self.oxides[ox]
                self.df[ph+'_mass'] = mass * self.df[ph]

        # Calculate the normalising value:
        for ph in phases_found:
            if ph is phases_found[0]:
                self.df['total_mass'] = self.df[ph+'_mass'].fillna(0.0)
            else:
                self.df['total_mass'] = self.df['total_mass'] + self.df[ph+'_mass'].fillna(0.0)

        # Normalise:
        for ph in phases_found:
            self.df[ph+'_mass'] = self.df[ph+'_mass']/self.df['total_mass']
        self.df['total_mass'] = self.df['total_mass']/self.df['total_mass']

        # Save the phases found:
        self.phases = phases_found
    
    def calculate_phase_wtpt(self):
        """
        Method for calculating the wtpt oxides of each phase.
        """
        for ph in self.phases:
            if ph + '_O' in self.df.columns:
                # First, convert FeO + O into FeO + Fe2O3:
                self.df[ph+'_Fe2O3'] = self.df[ph+'_O']
                self.df[ph+'_FeO'] = self.df[ph+'_FeO'] - 2*self.df[ph+'_O']
                self.df.drop(ph+'_O', axis=1, inplace=True)
        
            self.df[ph+'_total_mass'] = _np.zeros(self.df.shape[0])
            for ox in self.oxides:
                if ph + '_' + ox in self.df.columns:
                    self.df[ph + '_' + ox + '_wtpt'] = self.df[ph + '_' + ox] * self.oxides[ox]
                    self.df[ph + '_total_mass'] += self.df[ph + '_' + ox] * self.oxides[ox]
            
            for ox in self.oxides:
                if ph + '_' + ox in self.df.columns:
                    self.df[ph + '_' + ox + '_wtpt'] = self.df[ph + '_' + ox + '_wtpt'] / self.df[ph + '_total_mass'] * 100
            
            self.df.drop(ph + '_total_mass', axis=1, inplace=True)
    
    def calculate_phase_Mgn(self):
        """
        Calculate the Mg# of each phase from the mole fractions
        """
        for ph in self.phases:
            if ph + '_MgO' in self.df.columns and ph + '_FeO' in self.df.columns:
                self.df[ph + '_Mg#'] = self.df[ph + '_MgO'] / (self.df[ph + '_MgO'] + self.df[ph + '_FeO'])

    # def calculate_liquid_wtpt(self):
    #     """
    #     Method for converting mole fraction oxides to wtpt oxides.
    #     """
    #     # First, convert FeO + O into FeO + Fe2O3:
    #     self.df['liq_Fe2O3'] = self.df['liq_O']
    #     self.df['liq_FeO'] = self.df['liq_FeO'] - 2*self.df['liq_O']
    #     self.df.drop('liq_O', axis=1, inplace=True)

    #     self.df['liq_total_mass'] = _np.zeros(self.df.shape[0])
    #     for ox in self.oxides:
    #         if 'liq_' + ox in self.df.columns:
    #             self.df['liq_' + ox + '_wtpt'] = self.df['liq_' + ox] * self.oxides[ox]
    #             self.df['liq_total_mass'] += self.df['liq_' + ox] * self.oxides[ox]

    #     for ox in self.oxides:
    #         if 'liq_' + ox in self.df.columns:
    #             self.df['liq_' + ox + '_wtpt'] = self.df['liq_' + ox + '_wtpt'] / self.df['liq_total_mass'] * 100

    def make_f_grid(self, variable, deltaf=0.01, fill_value=[0.0, 0.0], use_edge_values=False):
        """
        Make a grid of a variable over a pressure- melt fraction grid. The scipy interpolate
        methods are used to construct a regularly spaced f grid at each supplied pressure. Where
        a variable becomes NaN, the change is extrapolated so that the variable will move
        continuously to zero (unless told not to do this).

        Parameters
        ----------
        variable : str
            The name of the variable in the phase diagram table
        deltaf : float, default: 0.01
            The melt fraction increment to use when generating the grid.
        fill_value : float or str, default: 0.0
            The fill_value to pass to scipy.interpolate.interpolate1d.
        use_edge_values : bool, default: False
            Instead of using zero or extrapolating for the edges of a defined variable- use the
            last finite value.

        Returns
        -------
        np.array
            3-dimensional array [0, :, :] is pressure, [1, :, :] is melt fraction, [2, :, :] is the
            value of the variable.
        """
        f = _np.arange(0.0, 1.0, deltaf)
        p = self.df['pressure'].unique()
        pp, ff = _np.meshgrid(p, f)
        grid = _np.zeros([3, _np.shape(pp)[0], _np.shape(pp)[1]])
        grid[0, :, :] = pp
        grid[1, :, :] = ff
        fill_value = _copy(fill_value)

        for i in range(_np.shape(pp)[1]):
            dfp = self.df[(self.df.pressure==pp[0, i])]
            x = dfp.liq_mass.copy()
            y = dfp[variable].copy()

            # Check that there's a liq=1 value, if not create one:
            if _np.isnan(x.iloc[len(x)-1]) == False and x.iloc[len(x)-1] < 1:
                ind = list(x.index)
                ind = ind + [ind[len(x)-1]+1]
                xnew = _np.zeros(_np.shape(x)[0]+1)
                xnew[:-1] = x
                xnew[-1] = 1.0
                x = _pd.Series(xnew, index=ind)
                ynew = _np.zeros(_np.shape(y)[0]+1)
                ynew[:-1] = y
                ynew[-1] = _np.nan
                y = _pd.Series(ynew, index=ind)

            # To ensure smoothness to zero, we will artificially overshoot F=0:
            first_finite_ind = x[(_np.isnan(x) == False)].index[0]
            dfdi0 = x.loc[first_finite_ind + 1] - x.loc[first_finite_ind]
            if dfdi0 < 0:
                raise Exception("Nonsense at"+ str(p[i]))
            # To accomodate phase diagrams where the first increments of melting are not present:
            dfdi_multiplier = 1
            while x.loc[first_finite_ind] - dfdi0 * dfdi_multiplier > 0:
                dfdi_multiplier += 1
            x.loc[first_finite_ind - 1] = x.loc[first_finite_ind] - dfdi0 * dfdi_multiplier

            # Now ensure that the variable can overshoot on both ends:
            y = y[(_np.isnan(x) == False)]

            if len(y[(_np.isnan(y) == False)]) > 1:
                # Check lower end:
                if _np.isnan(y.iloc[0]):
                    first_finite_ind = y[(_np.isnan(y) == False)].index[0]
                    dydx = (y.loc[first_finite_ind + 1] - y.loc[first_finite_ind])/(x.loc[first_finite_ind + 1] - x.loc[first_finite_ind])
                    y.loc[first_finite_ind - 1] = y.loc[first_finite_ind] + dydx*(x.loc[first_finite_ind - 1] - x.loc[first_finite_ind]) * dfdi_multiplier
                    if use_edge_values is True:
                        fill_value[0] = y.loc[first_finite_ind - 1]

                # Check upper end:
                if _np.isnan(y.iloc[len(y)-1]):
                    last_finite_ind = y[(_np.isnan(y) == False)].index[-1]
                    dydx = (y.loc[last_finite_ind] - y.loc[last_finite_ind - 1])/(x.loc[last_finite_ind] - x.loc[last_finite_ind - 1])
                    y.loc[last_finite_ind + 1] = y.loc[last_finite_ind] + dydx*(x.loc[last_finite_ind + 1] - x.loc[last_finite_ind]) #* dfdi_multiplier
                    if use_edge_values is True:
                        fill_value[1] = y.loc[last_finite_ind + 1]

                # The fill value must be a tuple for interp1d:
                if isinstance(fill_value, tuple):
                    fill_value_use = fill_value
                else:
                    fill_value_use = tuple(fill_value)


                x = x[(_np.isnan(x) == False) & (_np.isnan(y) == False)]

                if len(x) > 1:
                    y = y[(_np.isnan(x) == False) & (_np.isnan(y) == False)]

                    try:
                        pinterp = _interp1d(x, y, fill_value=fill_value_use, bounds_error=False, kind='quintic')
                    except:
                        try:
                            pinterp = _interp1d(x, y, fill_value=fill_value_use, bounds_error=False, kind='cubic')
                        except:
                            pinterp = _interp1d(x, y, fill_value=fill_value_use, bounds_error=False)

                    grid[2, :, i] = pinterp(f)

        return grid

class interpolate_grid(object):
    def __init__(self, grid, keep_positive=True):
        self.interp = _interp2d(grid[0,0], grid[1,:,0], grid[2], kind='cubic')
        self.keep_positive = keep_positive

    def __call__(self, x, y):
        value = self.interp(x, y)
        if self.keep_positive and value < 0:
            value = _np.array([0.0])
        return value

class phaseDiagram(object):
    """
    Reads in a phase diagram grid and sets up the spline interpolation functions.
    Access the interpolated functions by calling the class with the variable name
    and state as arguments.

    Parameters
    ----------
    grids : gridsThermocalc or gridsMelts instance
        The phaseDiagram grid, read in using either the gridsThermocalc or
        gridsMelts classes.
    variables : list or None, default: None
        The variables to add to the phaseDiagram. If None then all the columns
        in the grid will be added.
    extrapolate_to_zero : list, default: []
        Variables to extrapolate to zero. Phase mass fractions are always included
        in addition to the variables named in the list.
    """
    def __init__(self, grids, variables=None, extrapolate_to_zero=[]):
        self._interpolated_functions = {}

        self.P_range = [_np.nanmin(grids.df['pressure'].unique()),
                        _np.nanmax(grids.df['pressure'].unique())]

        for ph in grids.phases:
            extrapolate_to_zero.append(ph + '_mass')

        if variables is None:
            self.variables = grids.df.columns
        else:
            self.variables = variables

        for col in self.variables:
            if col not in ['pressure', 'X', 'liq', 'liq_mass']:
                if col not in extrapolate_to_zero:
                    grid = grids.make_f_grid(col, use_edge_values=True)
                else:
                    grid = grids.make_f_grid(col)
                self._interpolated_functions[col] = interpolate_grid(grid)

    def __call__(self, variable, state):
        return self._interpolated_functions[variable](state.P, state.F)[0]
    
    def plot_TxSection(self, parameter, F_range=[0,1], P_range=None, F_steps=100, P_steps=100,
                       mask_mineral_out=False):

        if P_range is None:
            P_range = self.P_range
        
        p = _np.linspace(P_range[0], P_range[1], P_steps)
        f = _np.linspace(F_range[0], F_range[1], F_steps)

        pp, ff, = _np.meshgrid(p, f)
        cc = _np.zeros(_np.shape(pp))

        f, a =_plt.subplots(figsize=(4.5,3.5))

        for i in range(_np.shape(pp)[0]):
            for j in range(_np.shape(pp)[1]):
                state = _pd.Series({'P': pp[i,j],'F': ff[i,j]})
                if mask_mineral_out is True and self(parameter[:3] + '_mass', state) < 1e-10:
                    cc[i,j] = _np.nan
                else:
                    cc[i,j] = self(parameter, state)

        cs = a.contourf(pp, ff, cc, levels=25, cmap=_plt.cm.Reds)

        cbar = f.colorbar(cs)
        cbar.set_label(parameter)

        a.set_xlabel('Pressure (GPa)')
        a.set_ylabel('Melt Fraction')


        return f, a



def load_phaseDiagram(name='thermocalc_klb1'):
    """
    Loads any of the phaseDiagrams included with pyMelt. At the moment this includes:
    - thermocalc_klb1
    - thermocalc_kg1
    - melts_klb1
    By default the function will return the THERMOCALC KLB1 phase diagram

    Parameters
    ----------
    name : str; default: 'thermocalc_klb1'

    Returns
    -------
    phaseDiagram
        The requested phaseDiagram object
    """

    available_models = ['thermocalc_klb1', 'thermocalc_kg1', 'melts_klb1']
    build_files = ['import_thermocalc_klb1.ipynb',
                   'import_thermocalc_kg1.ipynb',
                   'import_pmelts_klb1.ipynb']

    if name in available_models:
        pyMelt_path = os.path.dirname(os.path.realpath(__file__))
        try:
            f = open(pyMelt_path + '/phaseDiagrams/' + name +'.p', 'rb')

            phaseDiagram_object = pickle.load(f)

            f.close()

        except Exception:
            print("Building Phase Diagram Object...")

            original_directory = os.getcwd()

            env_dir = pyMelt_path + '/phaseDiagrams/build/'
            os.chdir(env_dir)

            filename = pyMelt_path + '/phaseDiagrams/build/' + build_files[available_models.index(name)]

            with open(filename) as fp:
                nb = load(fp)

            for cell in nb['cells']:
                if cell['cell_type'] == 'code':
                    source = ''.join(line for line in cell['source'] if not line.startswith('%'))
                    exec(source, globals(), locals())

            fp.close()

            os.chdir(original_directory)

            os.rename(pyMelt_path + '/phaseDiagrams/build/' + name + '.p', pyMelt_path + '/phaseDiagrams/' + name + '.p')

            f = open(pyMelt_path + '/phaseDiagrams/' + name +'.p')

            phaseDiagram_object = pickle.load(f)

            f.close()

        return phaseDiagram_object














