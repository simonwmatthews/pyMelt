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
    
    def calculateMineralProportions(self, method='invmel', species_objects=None, **kwargs):
        """
        Calculate the mineral proportions in the residue at each stage using one of the
        chemistry models. 

        NEED TO FINISH WRITING THIS DOCUMENTATION.
        """

        # Assemble the species_objects dictionary if not provided
        if species_objects is None:
            species_objects = {}
            for lith in self.mantle.names:
                method_recon = {}
                kwargs_recon = {}
                if isinstance(method, dict):
                    method_recon = method[lith]
                else:
                    method_recon = method
                for kw in kwargs:
                    if isinstance(kwargs[kw], dict):
                        kwargs_recon[kw] = kwargs[kw][lith]
                    else:
                        kwargs_recon[kw] = kwargs[kw]
                
                species_objects[lith] = self._create_species_objects({'Fake': _np.nan},
                                                                      method_recon,
                                                                      **kwargs_recon)
        
        for lith in self.mantle.names:

            # Prepare the variables for storage:
            mineralProps = None

            for i, row in self.lithologies[lith].iterrows():
                if row.F > 1e-15:
                    calc_return = species_objects[lith][0].mineralProportions(row)

                    if isinstance(calc_return, dict) is False:
                        raise InputError("This model does not support mineral proportion calculations")
                    
                    # Check if the table needs to be initialised:
                    if mineralProps is None:
                        mineralProps = _np.full([_np.shape(self.P)[0], len(calc_return)], _np.nan)
                        mineralNames = list(calc_return.keys())
                    
                    mineralProps[i, :] = list(calc_return.values())
            
            # Store this lithology's phase fractions:
            constructdf = _pd.DataFrame(mineralProps, columns=mineralNames)

            # Check if these phase fractions exist already:
            repeats = [value for value in mineralNames if value in self.composition[lith].columns]
            self.composition[lith].drop(repeats, inplace=True, axis=1)
            self.composition[lith] = _pd.concat([self.composition[lith], constructdf], axis=1)

            # Register the new entries:
            for mineral in mineralNames:
                self._composition_variable_type[lith][mineral] = 'mineralProportion'
    
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
                        label = ph + '_' + ox + '_wtpt'
                        if label in lith.phaseDiagram.variables:
                            labels.append(label)
                    nox = len(labels)
                    nsteps = len(self.P)
                    results = _np.zeros([nsteps, nox])
                    for i, row in self.lithologies[lithname].iterrows():
                        for j in range(nox):
                            label = labels[j]
                            results[i,j] = lith.phaseDiagram(label, row)


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
    
    def calculateTraceElements(self, cs, D=None, **kwargs):
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
        cs : nested dict or pandas.DataFrame
            The concentration of each element in each lithology prior to melting. If
            supplying a nested dict, then the first layer is the lithology names. If
            supplying a DataFrame then each lithology should have a row.
        D : dict, nested dict, pandas.DataFrame, or None, default: None
            The bulk partition coefficients for each element in each lithology. If a 
            single dict is supplied then the same partition coefficient will be used 
            for each lithology. If using mineral-melt partition coefficients then leave
            as None. The partition coefficient can be a float or a function which takes
            a state Series as input.
        
        Notes
        -----
        Explain the default mineral-melt partition coefficients etc.
        """
        for n in range(len(self.mantle.lithologies)):
            lith = self.mantle.lithologies[n]
            lithname = self.mantle.names[n]

            # Determine whether bulk partition coefficient are being used.
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
                    raise InputError("The input format for D was not recognised")
            
            if useBulk:
                # Do the calculation with bulk partition coefficients!
                if isinstance(cs, dict):
                    nel = len(cs[lithname])
                    elnames = list(cs[lithname].keys())
                    cslith = list(cs[lithname].values())
                elif isinstance(cs, _pd.DataFrame):
                    elnames = list(cs.columns)
                    cslith = list(cs.loc[lithname])
                    nel = len(elnames)
                
                nsteps = len(self.P)
                results = _np.zeros([nsteps, nel])
                for i in range(len(nel)):
                    el = elnames[i]
                    Del = Dlith[el]
                    for j, state in self.lithologies[lithname].iterrows():
                        if callable(Del):
                            Del_j = Del(state)
                        else:
                            Del_j = Del
                    
                    results[j, i] = 0.0

            else:
                # Do the calculation with mineral-melt partition coefficients!
                continue


        

    # def calculateChemistry_new(self, )


    def calculateChemistry(self, elements=None, species_objects=None, method='default', **kwargs):
        """
        DEPRECATED ROUTINE. NEEDS TO BE REMOVED.

        Calculate the composition of the melt according to default (or user defined) chemical
        models. The method may be run with default options if the mantle consists of only one
        lithology. Otherwise the parameters for each lithology must be specified, or pre-defined
        species objects must be provided.

        If elements is not set and the mantle is made of one lithology, the composition will be
        set to the depleted mantle of Workman & Hart (2005). The INVMEL forward model is used by
        default, except for Ba and Rb which are modelled using continuous_instantaneous. The
        mineral-specific partition coefficients for INVMEL are the constant values compiled by
        Gibson & Geist (2010). The bulk partition coefficients for Ba and Rb are from Workman &
        Hart (2005).

        Unless passing the species objects directly, the parameters used by the trace element
        model must be given also. If a constant value should be used for every lithology and
        every element (e.g., phi for continuous melting) it can be passed as a float. If a value
        varies for each element, but not each lithology, for example a mineral-specific D, a
        dict can be passed with each element name as a key. If the parameter varies for each
        element and lithology, it can be supplied as a nested-dictionary, much like the
        `elements` parameter. See Notes for a list of parameters for the default models.

        Parameters
        ----------
        elements : dict of dicts of floats (optional)
            A dictionary with each lithology as a key, the values are themselves dictionaries
            with species names as keys, and their concentrations as the values. The elements must
            be among the default element list. Ignore if supplying the species objects directly in
            the `species_objects` argument.
        species_objects : dict or None, default: None
            Each lithology name is a key of the dictionary, the values are a list of
            pyMelt.chemical_classes.species objects. If None, the species objects will be generated
            from the `elements` argument.
        method : string or dict or dict of dicts, default: 'default'
            One of 'default', 'batch', 'continuous_accumulated', 'continuous_instantaneous',
            'invmel', 'phase_diagram_trace', 'phase_diagram_major'. If using different models for
            different elements, specify them as a dictionary. This can be nested within another
            dictionary if you wish to use different combinations for each lithology.

        Notes
        -----

        The 'batch' melting routine uses:
         - D, the bulk partition coefficient

        The 'continuous' melting routine uses:
         - D, the bulk partition coefficient
         - phi, the porosity during melting (as a fraction, not percent).

        The 'invmel' melting routine uses:
         - olv_D, the olivine-melt partition coefficient
         - cpx_D, the clinopyroxene-melt partition coefficient
         - opx_D, the orthopyroxene-melt partition coefficient
         - spn_D, the spinel-melt partition coefficient
         - grt_D, the garnet-melt partition coefficient
         - plg_D, the plagioclase-melt partition coefficient
         - MineralProportions, the mass fraction of each mineral present (prior to melting) in the
           spinel, plagioclase, and garnet field.
         - cpxExhaustion, the melt fraction at which cpx (and grt/plg/spn) are exhausted.
         - garnetInCoeffs, coefficients controlling the P and T of the garnet in reaction
         - spinelOutCoeffs, coefficients controlling the P and T of the spinel out reaction
         - plagioclaseInInterval, The plagioclase in interval (in km).

        The 'phase_diagram_trace' melting routine uses:
         - FILL IN THE DOCUMENTATION!

        The 'phase_diagram_major' melting routine uses:
         - FILL IN THE DOCUMENTATION!
        """
        _warnings.warn("The fractional melting routine for phase diagrams is not implemented correctly.")

        # Check if using defaults, and assemble args if so:
        if method == 'default':
            default_kwargs = {'olv_D': _chemistry.data.olv_D,
                              'cpx_D': _chemistry.data.cpx_D,
                              'opx_D': _chemistry.data.opx_D,
                              'spn_D': _chemistry.data.spn_D,
                              'grt_D': _chemistry.data.grt_D,
                              'plg_D': _chemistry.data.plg_D,
                              'D': _chemistry.data.workman05_D}

            for argname in default_kwargs:
                if argname not in kwargs:
                    kwargs[argname] = default_kwargs[argname]

            method = _chemistry.default_methods

        # Check if all elements are provided for each lithology
        if elements is not None and len(elements) > 1:
            lithologies = list(elements.keys())
            if any(set(elements[lith].keys()) != set(elements[lithologies[0]].keys())
                   for lith in lithologies):
                elements = self._tidy_chemistry_inputs(elements)

        if elements is not None and any([el in ['P', 'T', 'F']
                                         for el in elements[list(elements.keys())[0]]]):
            raise InputError("Please rename elements so that none of P, T, or F appear.")

        # Assemble the species_objects dictionary if not provided
        if species_objects is None:
            if elements is None and self.mantle.number_lithologies == 1:
                print("Lithology composition is set to the depleted mantle of Workman & Hart "
                      "(2005).")
                elements = {self.mantle.names[0]: asdict(_chemistry.data.workman05_dmm())}
            elif elements is None:
                raise InputError("Either species_objects or elements must be provided.")

            species_objects = {}
            for lith in elements:
                method_recon = {}
                kwargs_recon = {}
                if (isinstance(method, dict)
                        and any(item in method.keys() for item in elements.keys())):
                    method_recon = method[lith]
                else:
                    method_recon = method
                for kw in kwargs:
                    if (isinstance(kwargs[kw], dict)
                            and any(item in kwargs[kw].keys() for item in elements.keys())):
                        kwargs_recon[kw] = kwargs[kw][lith]
                    else:
                        kwargs_recon[kw] = kwargs[kw]
                species_objects[lith] = self._create_species_objects(elements[lith],
                                                                     method_recon,
                                                                     **kwargs_recon)
                

        # Register the variables:
        for lith in species_objects:
            for species in species_objects[lith]:
                self._composition_variable_type[lith][species.name] = species.calculation_type

        # Check that the lithology names are correct
        for lith in species_objects:
            if lith not in self.mantle.names:
                raise InputError("The lithology specified (" + lith + ") was not found in "
                                 "the mantle from which this melting column was constructed.")

        # Iterate through calculations for each lithology:
        for lith in species_objects:

            # Prepare the storage for liquid compositions
            species_names = []
            for j in range(len(species_objects[lith])):
                species_names.append(species_objects[lith][j].name)
            results = _np.zeros([_np.shape(self.P)[0], len(species_objects[lith])])

            # Prepare the storage for solid + liquid compositions
            solidcomp = None
            mineral_element_labels = []

            for i, row in self.lithologies[lith].iterrows():
                for j in range(len(species_objects[lith])):
                    if row.F > 1e-15:
                        calc_return = species_objects[lith][j].composition(row)

                        # Check if method returns mineral compositions too:
                        if isinstance(calc_return, dict):
                            results[i, j] = calc_return['liq']
                            calc_return.pop('liq')
                            
                            # Check if the solid + liquid comp table needs to be initialised:
                            if solidcomp is None:
                                solidcomp = _np.full([_np.shape(self.P)[0], 
                                                       len(calc_return) * len(species_names)],
                                                       _np.NaN)
                                mineral_names = list(calc_return.keys())
                                for k in range(len(species_names)):
                                    mineral_element_labels += [mn + '_' + species_names[k] for mn in mineral_names]
                            
                            solidcomp[i, j*len(mineral_names) : (j+1)*len(mineral_names)] = list(calc_return.values())

                        # If the method just returns the liquid composition:
                        else:
                            results[i, j] = calc_return
                    else:
                        results[i, j] = _np.nan
                
            # Need to store this lithology's solid compositions, if generated:
            if solidcomp is not None:
                constructdf = _pd.DataFrame(solidcomp, columns=mineral_element_labels)

                # Register the variables:
                for label in mineral_element_labels:
                    if label.split('_')[0] != 'liq':
                        self._composition_variable_type[lith][label] = 'solidComposition'

                # Check if the element exists already:
                repeats = [value for value in mineral_element_labels if value in self.composition[lith].columns]
                self.composition[lith].drop(repeats, inplace=True, axis=1)
                self.composition[lith] = _pd.concat([self.composition[lith], constructdf], axis=1)


            # Store the liquid compositions
            constructdf = _pd.DataFrame(results, columns=species_names)

            # Checks if the element exists already, and drops in original if so:
            repeats = [value for value in species_names if value in self.composition[lith].columns]
            self.composition[lith].drop(repeats, inplace=True, axis=1)

            self.composition[lith] = _pd.concat([self.composition[lith], constructdf], axis=1)

    
    def calculateStableIsotopes(self, species, fractionationFactors, isotopeRatioLabel, 
                                bulk=0.0, fractionalExtraction=False, porosity=0.0, 
                                **kwargs):
        """
        Write some documentation here...

        Parameters
        ----------
        species: str 
            The species to calculate the stable isotope fractionation for. E.g., MgO.
            Must correspond to a species that has been calculated in the liquid and
            solid already.
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
        """

        _warnings.warn("Isotope ratios for solid phases where their phase fraction "
                       "goes < 0.01 are masked as a temporary fix to problematic "
                       "imports.")
        
        # Prepare the columns for results, and check the fractionationFactors input is correct
        if isinstance(fractionationFactors, dict):
            phases = list(fractionationFactors.keys())
            colnames = [isotopeRatioLabel]
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
            if species not in self.composition[lith].columns:
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
            self._composition_variable_type[lith][isotopeRatioLabel] = 'liquidIsotopeRatio' + suffix + '_' + species
            for ph in phases:
                self._composition_variable_type[lith][ph + '_' + isotopeRatioLabel] = 'solidIsotopeRatio_' + species

        for lith in self.mantle.names:

            results = _np.full([_np.shape(self.P)[0], len(phases) + 1], _np.nan)
            
            for i, row in self.composition[lith].iterrows():

                if row['F'] > 1e-15:

                    cliq = row[species]
                    xliq = row['F']

                    # Assemble arrays for the summations:
                    x = _np.zeros(len(fractionationFactors))
                    c = _np.zeros(len(fractionationFactors))
                    a = _np.zeros(len(fractionationFactors))
                    
                    for n in range(len(phases)):
                        x[n] = row[phases[n]] * (1.0-xliq)
                        c[n] = row[phases[n] + '_' + species]
                        if callable(fractionationFactors[phases[n]]):
                            a[n] = fractionationFactors[phases[n]](row)
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

                        if F_prev > 1e-15:

                            bulk_a = (_np.sum(x * c * a) + porosity * cliq) / (_np.sum(x * c) + porosity * cliq)

                            # We might expect this block should be calculated in the previous step, and so we should
                            # be able to use it, except when I tried it completely broke the code. I think this is
                            # because what we need in this step is the final composition of melt extracted, rather
                            # than the aggregate, which is what is calculated in the previous step.
                            cs_prev_denom = 0.0
                            for ph in phases:
                                cs_prev_denom += row_prev[ph + '_' + species] * row_prev[ph]
                            cs_prev_numer = cs_prev_denom * (bulk[lith] / 1e3 + 1)
                            D_prev_denom = cs_prev_denom / row_prev[species]
                            D_prev_numer = D_prev_denom * bulk_a
                            
                            cs_denom = 0.0
                            for ph in phases:
                                cs_denom += row[ph + '_' + species] * row[ph]
                            D_denom = cs_denom / row[species]

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
