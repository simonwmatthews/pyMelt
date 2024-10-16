"""
===============
Melting Columns
===============

the meltingcolumn_classes module provides the melting column classes. At the moment it consists of
a single melting column- a simple 1D melting column.

"""

import numpy as _np
from pyMelt.core import InputError
import matplotlib.pyplot as _plt
import pandas as _pd
import pyMelt.chemistry as _chemistry
from dataclasses import asdict
from copy import copy as _copy
import warnings as _warnings


class meltingColumn():
    """
    Class for storing the results of a 1D multi-lithology melting model.

    Parameters
    ----------
    calculation_results : pandas.DataFrame
        Dataframe with columns 'P' for Pressure in GPa, 'T' for Temperature in
        degrees C, Remaining columns for melt fraction from each lithology.
    mantle : pyMelt.Mantle
        The mantle object used to generate the melting column.
    Tp : float
        The potential temperature used to generate the melting column, if applicable.

    Attributes
    ----------
    calculation_results : pandas.DataFrame
        The stored raw calculation results.
    mantle : pyMelt.Mantle
        The mantle object used to generate the melting column.
    Tp : float
        The potential temperature used to generate the melting column, if applicable.
    P : pandas.Series
        The pressure steps in GPa
    T : pandas.Series
        The temperature steps in degC.

    """

    def __init__(self, calculation_results, mantle, Tp=None):
        self.calculation_results = calculation_results
        self.mantle = mantle
        self.Tp = Tp
        self.P = calculation_results.P
        self.T = calculation_results['T']
        self.chemistry_output = None


        # New composition registry structure
        # Keeps track of what kind of variable each column in the composition tables
        # is, i.e., concentration, mineral composition, stable isotope ratio...
        self._composition_variable_type = {}
        for lith in self.mantle.names:
            self._composition_variable_type[lith] = {'P': 'metadata', 
                                                     'T': 'metadata',
                                                     'F': 'metadata'}

        # New structure for storing composition information (solid + liquid)
        self.composition = {}

        self.lithologies = {}
        for i in range(self.mantle.number_lithologies):
            df = _pd.DataFrame()
            df['P'] = self.P
            df['T'] = self.T
            df['F'] = self.calculation_results[self.mantle.names[i]]
            self.lithologies[self.mantle.names[i]] = df
            self.composition[self.mantle.names[i]] = _copy(df) # These contributions to the composition table are already registered above
            
        self.F = _np.zeros(_np.shape(self.P)[0])

        for i in range(self.mantle.number_lithologies):
            self.F = (self.F + self.mantle.proportions[i]
                      * self.calculation_results[self.mantle.names[i]])

    def plot(self):
        """
        Generates a plot showing the thermal gradient and melt fractions of each lithology.

        Returns
        -------
        (matplotlib.figure, matplotlib.axes)
            The generated figure and axes objects.
        """
        f, a = _plt.subplots(1, 2, sharey='row', dpi=100)

        lith = self.mantle.names

        for i in range(_np.shape(lith)[0]):
            a[1].plot(self.lithologies[lith[i]].F, self.P, label=lith[i])

        a[0].plot(self.T, self.P, label='Thermal Gradient', c='k')
        a[1].plot(self.F, self.P, label='Total', c='k', ls='--')

        P = _np.linspace(_np.min(self.P), _np.max(self.P), 101)
        for i in range(self.mantle.number_lithologies):
            T = self.mantle.lithologies[i].TSolidus(P)
            a[0].plot(T, P, label=self.mantle.names[i] + ' solidus')

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
    
    def calculateMineralProportions(self):
        """
        Calculate the mineral proportions in the residue at each calculation step using the
        phase diagrams attached to each lithology. The results are stored in each lithology's
        table in meltingColumn.composition.
        """

        for n in range(len(self.mantle.names)):
            lith = self.mantle.lithologies[n]
            lithname = self.mantle.names[n]

            if lith.phaseDiagram is None:
                raise InputError("A phase diagram must be attached to the lithology in order"
                                 " to calculate mineral proportions.")
            
            nmin = len(lith.phaseDiagram.minerals)
            nsteps = len(self.P)
            minprops = _np.zeros([nsteps, nmin])

            for i, state in self.lithologies[lithname].iterrows():
                for j in range(nmin):  
                    minprops[i,j] = lith.phaseDiagram(lith.phaseDiagram.minerals[j] + "_mass", state)
            
            # Store this lithology's phase fractions:
            constructdf = _pd.DataFrame(minprops, columns=lith.phaseDiagram.minerals)

            # Check if these phase fractions exist already:
            repeats = [value for value in lith.phaseDiagram.minerals if value in self.composition[lithname].columns]
            self.composition[lithname].drop(repeats, inplace=True, axis=1)
            self.composition[lithname] = _pd.concat([self.composition[lithname], constructdf], axis=1)

            # Register the new entries:
            for mineral in lith.phaseDiagram.minerals:
                self._composition_variable_type[lithname][mineral] = 'mineralProportion'
            
    
    def calculateMajorOxides(self):
        """
        Calculate the major element oxide concentrations of the melts generated from each
        lithology. The concentrations are looked up from the phase diagrams attached to
        each lithology. The results are stored in column.composition.
        """

        for n in range(len(self.mantle.lithologies)):
            lith = self.mantle.lithologies[n]
            lithname = self.mantle.names[n]
            if lith.phaseDiagram is not None:
                for ph in ['liq'] + lith.phaseDiagram.minerals:
                    labels = []
                    # Check which oxides are in this phase:
                    for ox in lith.phaseDiagram.oxides:
                        label = ph + '_' + ox #+ '_wtpt'
                        if label + '_wtpt' in lith.phaseDiagram.variables:
                            labels.append(label)
                    nox = len(labels)
                    nsteps = len(self.P)
                    results = _np.zeros([nsteps, nox])
                    for i, row in self.lithologies[lithname].iterrows():
                        for j in range(nox):
                            label = labels[j]
                            results[i,j] = lith.phaseDiagram(label + '_wtpt', row)


                    constructdf = _pd.DataFrame(results, columns=labels)

                    # Register the variables:
                    for label in labels:
                        if label.split('_')[0] != 'liq':
                            self._composition_variable_type[lithname][label] = 'solidComposition'
                        else:
                            self._composition_variable_type[lithname][label] = "liquidConcentrationInstantaneous"
                    
                    # Check if the element exists already:
                    repeats = [value for value in labels if value in self.composition[lithname].columns]
                    self.composition[lithname].drop(repeats, inplace=True, axis=1)
                    self.composition[lithname] = _pd.concat([self.composition[lithname], constructdf], axis=1)
    
    def calculateTraceElements(self, c0, D=None, porosity=0.0, calcMineralCompositions=True, **kwargs):
        """
        Calculate the trace element composition of liquid (and solid) using partition
        coefficients. These partition coefficients may represent the bulk partitioning
        (by setting D) or mineral-melt partitioning by supplying arguments formatted as
        D_min. If supplying mineral-melt partition coefficients then the lithologies must
        have a phaseDiagram attached to them.

        Default mineral-melt partition coefficients are automatically used if none
        are specified. See Notes.

        Parameters
        ----------
        c0 : nested dict or pandas.DataFrame
            The concentration of each element in each lithology prior to melting. If
            supplying a nested dict, then the first layer is the lithology names. If
            supplying a DataFrame then each lithology should have a row.
        D : dict, nested dict, pandas.DataFrame, or None, default: None
            The bulk partition coefficients for each element in each lithology. If a 
            single dict is supplied then the same partition coefficient will be used 
            for each lithology. If using mineral-melt partition coefficients then leave
            as None. The partition coefficient can be a float or a function which takes
            a state Series as input.
        porosity: float or dict, default: 0.0
            The residual porosity during melting, as a fraction. Different values can
            be set for each lithology by using a dict, with the lithology name as the key.
        calcMineralCompositions: bool, default: True
            If mineral-melt partition coefficients are supplied then calculate and store
            the trace element concentrations in each mineral too. You might want to disable
            this if you need the calculation to run as fast as possible and you are only
            interested in the liquid composition.
        
        Notes
        -----
        Explain the default mineral-melt partition coefficients etc.
        """

        # Set up default partition coefficients if not supplied:
        if "D_olv" not in kwargs:
            kwargs['D_olv'] = asdict(_chemistry.data.olv_D())
        if "D_cpx" not in kwargs:
            kwargs['D_cpx'] = asdict(_chemistry.data.cpx_D())
        if "D_opx" not in kwargs:
            kwargs['D_opx'] = asdict(_chemistry.data.opx_D())
        if "D_grt" not in kwargs:
            kwargs['D_grt'] = asdict(_chemistry.data.grt_D())
        if "D_spn" not in kwargs:
            kwargs['D_spn'] = asdict(_chemistry.data.spn_D())
        if "D_plg" not in kwargs:
            kwargs['D_plg'] = asdict(_chemistry.data.plg_D())

        # Handle one lithology at a time to allow for different number of calculation steps
        for n in range(len(self.mantle.lithologies)):
            lith = self.mantle.lithologies[n]
            lithname = self.mantle.names[n]

            # Get the porosity
            if isinstance(porosity, dict):
                lithporosity = porosity[lithname]
            elif isinstance(porosity, float):
                lithporosity = porosity
            else:
                raise InputError("The input format for porosity was not recognised.")

            # Determine whether bulk partition coefficients are being used.
            useBulk = False
            if D is not None:
                if isinstance(D, dict):
                    if isinstance(D[list(D.values())[0]], dict):
                        # Nested dictionary
                        if lithname in D:
                            # Nested dictionary that has this lithology!
                            Dlith = D[lithname]
                            useBulk = True
                            continue
                    else:
                        # Non-nested dict, so use same D for all lith
                        Dlith = D
                        useBulk = True
                elif isinstance(D, _pd.DataFrame):
                    if lithname in D.index:
                        Dlith = D.loc[lithname]
                        useBulk = True
                else:
                    raise InputError("The input format for D was not recognised.")
            
            # Prepare calculation inputs
            if isinstance(c0, dict):
                nel = len(c0[lithname])
                elnames = list(c0[lithname].keys())
                cs = _np.array(list(c0[lithname].values()))
            elif isinstance(c0, _pd.DataFrame):
                elnames = list(c0.columns)
                cs = _np.array(list(c0.loc[lithname]))
                nel = len(elnames)
            
            
            # Setup output array
            nsteps = len(self.P)
            results = _np.zeros([nsteps, nel])

            # Retrieve bulk D and P for each step
            prevState = None

            # Create an array to store mineral-melt partition coefficients in case they are needed
            # to calculate mineral compositions
            mineralDs = _np.zeros([nsteps, nel, len(lith.phaseDiagram.minerals)])

            firstcalc = True
            for i, state in self.composition[lithname].iterrows():
                if prevState is None and state['F'] > 1e-15:
                    prevState = _copy(state)
                    prevState['F'] = 0.0
                # Don't calculate using the first step, as we need a change in min props to calculate P
                if state['F'] > 1e-15:
                # if prevState is not None and state['F'] > 1e-15:

                    # Assemble bulk partition coefficients
                    bulkD = _np.zeros([nel])
                    bulkP = _np.zeros([nel])
                    if useBulk:
                    # Do the calculation with bulk partition coefficients!
                        for j in range(len(nel)):
                            el = elnames[j]
                            Del = Dlith[el]
                            if callable(Del):
                                bulkD[j] = (Del(state, lith, **kwargs) + lithporosity) / (1 + lithporosity)
                                bulkP[j] = bulkD[i,j] # If D is constant then P is the same
                            else:
                                bulkD[j] = (Del + lithporosity) / (1 + lithporosity)
                                bulkP[j] = Del # If D is constant then P is the same

                    else:
                    # Calculate bulk partition coefficients from mineral proportions

                        for j in range(nel):
                            el = elnames[j]
                            for k in range(len(lith.phaseDiagram.minerals)):
                                min = lith.phaseDiagram.minerals[k]
                                if "D_" + min in kwargs:
                                    Del = kwargs["D_"+min][el]
                                    if callable(Del):
                                        stepD = Del(state, lith, **kwargs)
                                    else:
                                        stepD = Del
                                    mineralDs[i, j, k] = stepD
                                    bulkD[j] += state[min] * stepD
                                    if firstcalc:
                                        bulkP[j] = bulkD[j]
                                    else:
                                        bulkP[j] += (prevState[min] * (1-prevState['F']) - state[min] * (1 - state['F']) ) * stepD / (1- state['F'])
                            bulkD[j] = (bulkD[j] + lithporosity) / (1 + lithporosity)
                    
                    # Integrate Shaw equations to find new cs and cl
                    cl = cs * (1 - prevState['F']) / (bulkD - bulkP * prevState['F'])
                    k1 = (cs - cl) / (1 - prevState['F'])
                    
                    F = prevState['F'] + (state['F'] - prevState['F']) / 2
                    cs = cs + k1 * (state['F'] - prevState['F']) / 2
                    cl = cs * (1 - F) / (bulkD - bulkP * F)
                    k2 = (cs - cl) / (1 - F)

                    # Same F as for k2
                    cs = cs + k2 * (state['F'] - prevState['F']) / 2
                    cl = cs * (1 - F) / (bulkD - bulkP * F)
                    k3 = (cs - cl) / (1 - F)

                    F = state['F']
                    cs = cs + k3 * (state['F'] - prevState['F'])
                    cl = cs * (1 - F) / (bulkD - bulkP * F)
                    k4 = (cs - cl) / (1 - F)

                    cs = cs + (1 / 6) * (state['F'] - prevState['F']) * (k1 + 2*k2 + 2*k3 + k4)
                    cl = cs * (1 - F) / (bulkD - bulkP * F)

                    # Store results
                    results[i, :] = cl
                    firstcalc = False
                else:
                    results[i, :] = _np.nan

                # Store this state for use in the next step of the calculation
                prevState = state

            # Store this in the composition dataframe
            colnames = []
            for i in range(nel):
                colnames.append("liq_" + elnames[i])

            constructdf = _pd.DataFrame(results, columns=colnames)
            # Check if the element exists already:
            repeats = [value for value in colnames if value in self.composition[lithname].columns]
            self.composition[lithname].drop(repeats, inplace=True, axis=1)
            self.composition[lithname] = _pd.concat([self.composition[lithname], constructdf], axis=1)

            # Register the new entries:
            for colname in colnames:
                self._composition_variable_type[lithname][colname] = 'liquidConcentrationInstantaneous'
            

            # Calculate mineral compositions if required
            if calcMineralCompositions and useBulk is False:
                for k in range(len(lith.phaseDiagram.minerals)):
                    min = lith.phaseDiagram.minerals[k]
                    mincomp = results * mineralDs[:, :, k]
                    colnames = [min + "_" + elname for elname in elnames]
                    constructdf = _pd.DataFrame(mincomp, columns=colnames)

                    # Check if the element exists already:
                    repeats = [value for value in colnames if value in self.composition[lithname].columns]
                    self.composition[lithname].drop(repeats, inplace=True, axis=1)
                    self.composition[lithname] = _pd.concat([self.composition[lithname], constructdf], axis=1)

                    # Register the new entries:
                    for colname in colnames:
                        self._composition_variable_type[lithname][colname] = 'mineralConcentration'


    
    def calculateStableIsotopes(self, species, fractionationFactors, isotopeRatioLabel, 
                                bulk=0.0, fractionalExtraction=False, porosity=0.0, 
                                speciesMassConversion=[],
                                **kwargs):
        """
        Write some documentation here...

        Parameters
        ----------
        species: str , list
            The species to calculate the stable isotope fractionation for. E.g., MgO.
            Must correspond to a species that has been calculated in the liquid and
            solid already. If multiple species exist for the same isotope system (e.g.,
            FeO and Fe2O3), then provide this as a list, and specify a species_mass_conversion.
        fractionationFactors : dict
            The mineral-liquid fractionation factors (1000 ln beta) for each mineral 
            in the calculation. The mineral name should be given as the key. A number 
            or a function may be given as the value. 
        isotopeRatioLabel : str
            The label to be applied to the results, e.g., 'd57Fe'.
        bulk : float or dict, default: 0
            The bulk isotope ratio, in the units of the calculation. If the calculation
            is simulating fractional melt extraction this number corresponds to the
            bulk value before melting. Use a dictionary to provide different bulk
            compositions for each lithology (with the lithology name as the key).
        fractionalExtraction : bool, default: False
            Controls wether calculation assumes batch or fractional melting for the
            purposes of the isotope fractionation calculation.
        porosity : float, default: 0.0
            If doing a fractional melting calculation, phi allows some residual melt to be
            retained at each step (i.e., continuous melting). Perfect fractional melting
            is assumed by default. If modelling a trace element system then this should
            be set to the same value used during the calculation of trace element
            concentrations.
        speciesMassConversion : list, default: []
            If using multiple species (e.g., FeO and Fe2O3), then specify the constant with
            which to divide each by in order to conserve mass when the two variables are
            combined. If this is FeO and Fe2O3, then this will be set automatically.
        """

        _warnings.warn("Isotope ratios for solid phases where their phase fraction "
                       "goes < 0.01 are masked as a temporary fix to problematic "
                       "imports.")
        
        # If the species has multiple forms, e.g., FeO + Fe2O3
        # There may be a more elegant way of doing this, but this should work for now.

        if isinstance(species, list):
            if len(speciesMassConversion) == 0:
                if (species[0] == 'FeO' and species[1] == 'Fe2O3'):
                    speciesMassConversion = [1.0, 1.111]
                elif (species[0] == 'Fe2O3' and species[1] == 'FeO'):
                    speciesMassConversion = [1.111, 1.0]
                else:
                    raise InputError("Multiple species declared, but no mass conversion specified.")

            phases = list(fractionationFactors.keys()) + ['liq']

            specieslabel = species[0] + 'T'

            for lith in self.mantle.names:
                for ph in phases:
                    totspeciesconc = _np.zeros(len(self.composition[lith]))
                    for n in range(len(species)):
                        sp = species[n]
                        totspeciesconc += self.composition[lith][ph + "_" + sp] / speciesMassConversion[n]
                    self.composition[lith][ph + '_' + specieslabel] = totspeciesconc
                    self._composition_variable_type[lith][ph + '_' + specieslabel] = 'liquidConcentrationInstantaneous'
            species = specieslabel
        
        # Prepare the columns for results, and check the fractionationFactors input is correct
        if isinstance(fractionationFactors, dict):
            phases = list(fractionationFactors.keys())
            colnames = ['liq_' + isotopeRatioLabel]
            for ph in phases:
                colnames.append(ph + '_' + isotopeRatioLabel)
        else:
            raise InputError("fractionationFactors must be a dict. If you want to use a single "
                            "fractionation factor for solid-liquid fractionation then specify "
                            "the same value for each mineral.")
        
        # Check the bulk composition input
        if isinstance(bulk, float): 
            if len(self.mantle.names) > 1:
                _warnings.warn("A single bulk isotope ratio is being applied to every lithology. "
                               "Unless there is no isotopic heterogeneity this means the calculation "
                               "will only be indicative of general behaviour.")
            bulkval = bulk
            bulk = {}
            for lith in self.mantle.names:
                bulk[lith] = bulkval
        else:
            for bn in bulk:
                if bn not in self.mantle.names:
                    raise InputError("The lithology {} was not recognised.".format(bn))
        
        # Check the species exists for each of the phases:
        for lith in self.mantle.names:
            if 'liq_' + species not in self.composition[lith].columns:
                raise InputError("{0} was not found in {1} for {2}. The composition of each "
                                    "phase must have already been calculated.".format(species, 'liquid', lith))
            for ph in phases:
                if ph + '_' + species not in self.composition[lith].columns:
                    raise InputError("{0} was not found in {1} for {2}. The composition of each "
                                    "phase must have already been calculated.".format(species, ph, lith))
        
        
        
        # Register the variables:
        if fractionalExtraction is True:
            suffix = 'Instantaneous'
        else:
            suffix = 'Aggregated'
        for lith in self.mantle.names:
            self._composition_variable_type[lith]['liq_' + isotopeRatioLabel] = 'liquidIsotopeRatio' + suffix + '_' + species
            for ph in phases:
                self._composition_variable_type[lith][ph + '_' + isotopeRatioLabel] = 'solidIsotopeRatio_' + species

        for n in range(len(self.mantle.names)):
            lithobj = self.mantle.lithologies[n]
            lith = self.mantle.names[n]

            results = _np.full([_np.shape(self.P)[0], len(phases) + 1], _np.nan)
            
            for i, row in self.composition[lith].iterrows():

                if row['F'] > 1e-15:

                    cliq = row['liq_' + species]
                    xliq = row['F']

                    # Assemble arrays for the summations:
                    x = _np.zeros(len(fractionationFactors))
                    c = _np.zeros(len(fractionationFactors))
                    a = _np.zeros(len(fractionationFactors))
                    
                    for n in range(len(phases)):
                        x[n] = row[phases[n]] * (1.0-xliq)
                        c[n] = row[phases[n] + '_' + species]
                        if callable(fractionationFactors[phases[n]]):
                            a[n] = fractionationFactors[phases[n]](row, lithobj.phaseDiagram)
                        else:
                            a[n] = fractionationFactors[phases[n]]
                    a = _np.exp(a/1000)
                                            
                
                    if fractionalExtraction is False:

                        cbulk = _np.sum(c * x) + cliq * xliq

                        # delta_melt = (
                        #     (cbulk * (bulk[lith]/1e3 + 1) 
                        #     / (np.sum(x*c*a) + xliq * cliq)
                        #     - 1) *1e3
                        # )

                        delta_melt = (
                                (cbulk * bulk[lith] - 1e3 * (_np.sum(x*c*a) - _np.sum(x*c)))
                                / (_np.sum(x*c*a) + xliq * cliq)
                            )
                        
                        results[i, 0] = delta_melt
                        results[i, 1:] = a * (delta_melt + 1e3) - 1e3

                        for j in range(len(x)):
                            if x[j] < 1e-2:
                                results[i, j+1] = _np.nan
                    
                    else:
                    
                        row_prev = self.composition[lith].iloc[i-1]
                        F_prev = row_prev['F']

                        if F_prev > 1e-15 and row_prev['liq_'+species] > 0:

                            bulk_a = (_np.sum(x * c * a) + porosity * cliq) / (_np.sum(x * c) + porosity * cliq)

                            # We might expect this block should be calculated in the previous step, and so we should
                            # be able to use it, except when I tried it completely broke the code. I think this is
                            # because what we need in this step is the final composition of melt extracted, rather
                            # than the aggregate, which is what is calculated in the previous step.
                            cs_prev_denom = 0.0
                            for ph in phases:
                                cs_prev_denom += row_prev[ph + '_' + species] * row_prev[ph]
                            cs_prev_numer = cs_prev_denom * (bulk[lith] / 1e3 + 1)
                            D_prev_denom = cs_prev_denom / row_prev['liq_' + species]
                            D_prev_numer = D_prev_denom * bulk_a
                            
                            cs_denom = 0.0
                            for ph in phases:
                                cs_denom += row[ph + '_' + species] * row[ph]
                            D_denom = cs_denom / row['liq_' + species]

                            D_numer = D_denom * bulk_a

                            cl_prev_denom = cs_prev_denom / D_prev_denom 
                            cl_prev_numer = cs_prev_numer / D_prev_numer
                            # End of block of maybe repeated/redundant code.
                            
                            norm = 0.0
                            pbar_numer = _np.zeros(len(phases))
                            pbar_denom = _np.zeros(len(phases))
                            for j in range(len(phases)):
                                ph = phases[j]
                                norm += (row_prev[ph] * (1-F_prev) - row[ph] * (1-row['F']))
                                pbar_numer[j] = (row_prev[ph] * (1-F_prev) - row[ph] * (1-row['F'])) * row[ph+'_'+species] / cliq * a[j]
                                pbar_denom[j] = (row_prev[ph] * (1-F_prev) - row[ph] * (1-row['F'])) * row[ph+'_'+species] / cliq
                            Pbar_numer = _np.sum(pbar_numer) / norm
                            Pbar_denom = _np.sum(pbar_denom) / norm

                            # Calculate dcs/dX over integration range (Shaw eqns)
                            k1_numer = (cs_prev_numer - cl_prev_numer) / (1 - F_prev)
                            k1_denom = (cs_prev_denom - cl_prev_denom) / (1 - F_prev)


                            dF = (row['F'] - F_prev) / 2
                            k_cs = cs_prev_numer + k1_numer * dF
                            k_cl = k_cs * (1 - F_prev - dF) / (D_numer * (1 - F_prev) - Pbar_numer * (dF))
                            k2_numer = (k_cs - k_cl) / (1 - (F_prev + dF))
                            k_cs = cs_prev_denom + k1_denom * dF
                            k_cl = k_cs * (1 - F_prev - dF) / (D_denom * (1 - F_prev) - Pbar_denom * (dF))
                            k2_denom = (k_cs - k_cl) / (1 - (F_prev + dF))

                            dF = (row['F'] - F_prev) / 2
                            k_cs = cs_prev_numer + k2_numer * dF
                            k_cl = k_cs * (1 - F_prev - dF) / (D_numer * (1 - F_prev) - Pbar_numer * (dF))
                            k3_numer = (k_cs - k_cl) / (1 - (F_prev + dF))
                            k_cs = cs_prev_denom + k2_denom * dF
                            k_cl = k_cs * (1 - F_prev - dF) / (D_denom * (1 - F_prev) - Pbar_denom * (dF))
                            k3_denom = (k_cs - k_cl) / (1 - (F_prev + dF))

                            dF = (row['F'] - F_prev)
                            k_cs = cs_prev_numer + k3_numer * dF
                            k_cl = k_cs * (1 - F_prev - dF) / (D_numer * (1 - F_prev) - Pbar_numer * (dF))
                            k4_numer = (k_cs - k_cl) / (1 - (F_prev + dF))
                            k_cs = cs_prev_denom + k3_denom * dF
                            k_cl = k_cs * (1 - F_prev - dF) / (D_denom * (1 - F_prev) - Pbar_denom * (dF))
                            k4_denom = (k_cs - k_cl) / (1 - (F_prev + dF))
                            
                            
                            cs_numer = cs_prev_numer + (1 / 6) * (row['F'] - F_prev) * (k1_numer + 2*k2_numer + 2*k3_numer + k4_numer)
                            cl_numer = cs_numer * (1 - row['F']) / (D_numer * (1 - F_prev) - Pbar_numer * (row['F'] - F_prev))
                            cs_denom = cs_prev_denom + (1 / 6) * (row['F'] - F_prev) * (k1_denom + 2*k2_denom + 2*k3_denom + k4_denom)
                            cl_denom = cs_denom * (1 - row['F']) / (D_denom * (1 - F_prev) - Pbar_denom * (row['F'] - F_prev))


                            bulk[lith] = (cs_numer / cs_denom - 1) * 1e3
                            results[i, 0] = (cl_numer / cl_denom - 1) * 1e3
                            results[i, 1:] = a * (results[i,0] + 1e3) - 1e3

                            for j in range(len(x)):
                                if x[j] < 1e-2:
                                    results[i, j+1] = _np.nan
                        

            constructdf = _pd.DataFrame(results, columns=colnames)
            # Check if the element exists already:
            repeats = [value for value in colnames if value in self.composition[lith].columns]
            self.composition[lith].drop(repeats, inplace=True, axis=1)
            self.composition[lith] = _pd.concat([self.composition[lith], constructdf], axis=1)

    def _create_species_objects(self, elements, method, **kwargs):
        methods = {'batch': _chemistry.batchSpecies,
                   'continuous_instantaneous': _chemistry.continuousSpecies_instantaneous,
                   'continuous_accumulated': _chemistry.continuousSpecies_accumulated,
                   'invmel': _chemistry.invmelSpecies,
                   'phase_diagram_trace': _chemistry.phaseDiagramTraceSpecies,
                   'phase_diagram_major': _chemistry.phaseDiagramMajorSpecies}
        species_objects = []
        for el in elements:
            kwargs_recon = {}
            for arg in kwargs:
                if (isinstance(kwargs[arg], dict)
                        and any(item in kwargs[arg].keys() for item in elements.keys())):
                    if el in kwargs[arg].keys():
                        kwargs_recon[arg] = kwargs[arg][el]
                    else:
                        kwargs_recon[arg] = None
                else:
                    kwargs_recon[arg] = kwargs[arg]
            if isinstance(method, dict) and any(item in method.keys() for item in elements.keys()):
                species_objects.append(methods[method[el]](el, c0=elements[el], **kwargs_recon))
            else:
                species_objects.append(methods[method](el, c0=elements[el], **kwargs_recon))
        return species_objects

    def _tidy_chemistry_inputs(self, elements):
        lithologies = list(elements.keys())

        elements_list = list(elements[lithologies[0]].keys())
        for lith in lithologies:
            to_pop = []
            for i in range(len(elements_list)):
                if elements_list[i] not in elements[lith]:
                    to_pop.append(i)
            for i in range(len(to_pop)):
                elements_list.pop(to_pop[i] - i)

        reconstruct = {}
        for lith in lithologies:
            reconstruct[lith] = {el: elements[lith][el] for el in elements_list}

        return reconstruct
