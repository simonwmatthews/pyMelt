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
import pandas as pd
import pyMelt.chemistry
from copy import copy as _copy


class meltingColumn():
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
    calculation_results : pandas.DataFrame
        The stored raw calculation results.
    mantle : pyMelt.Mantle
        The mantle class used to generate the melting column.
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
        self._species_calc_type = None

        self.lithologies = {}
        for i in range(self.mantle.number_lithologies):
            df = pd.DataFrame()
            df['P'] = self.P
            df['T'] = self.T
            df['F'] = self.calculation_results[self.mantle.names[i]]
            self.lithologies[self.mantle.names[i]] = df

        self.F = np.zeros(np.shape(self.P)[0])

        for i in range(self.mantle.number_lithologies):
            self.F = (self.F + self.mantle.proportions[i]
                      * self.calculation_results[self.mantle.names[i]])

    def plot(self):
        """
        Generates a plot showing the thermal gradient and melt fractions ofeach lithology.

        Returns
        -------
        (matplotlib.figure, matplotlib.axes)
            The generated figure and axes objects.
        """
        f, a = plt.subplots(1, 2, sharey='row', dpi=100)

        lith = self.mantle.names

        for i in range(np.shape(lith)[0]):
            a[1].plot(self.lithologies[lith[i]].F, self.P, label=lith[i])

        a[0].plot(self.T, self.P, label='Thermal Gradient', c='k')
        a[1].plot(self.F, self.P, label='Total', c='k', ls='--')

        P = np.linspace(np.min(self.P), np.max(self.P), 101)
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

    def calculateChemistry(self, elements=None, species_objects=None, method='default', **kwargs):
        """
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
            'invmel'. If using different models for different elements, specify them as a
            dictionary. This can be nested within another dictionary if you wish to use different
            combinations for each lithology.

        Notes
        -----

        The 'batch' melting routine uses:
         - D, the bulk partition coefficient

        The 'continuous' melting routine uses:
         - D, the bulk partition coefficient
         - phi, the porosity during melting (as a fraction, not percent).

        The 'invmel' melting routine uses:
         - olv_D, the olivine-melt partition coefficient
         - cpx_D, the olivine-melt partition coefficient
         - opx_D, the olivine-melt partition coefficient
         - spn_D, the olivine-melt partition coefficient
         - grt_D, the olivine-melt partition coefficient
         - plg_D, the olivine-melt partition coefficient
         - MineralProportions, the mass fraction of each mineral present (prior to melting) in the
           spinel, plagioclase, and garnet field.
         - cpxExhaustion, the melt fraction at which cpx (and grt/plg/spn) are exhausted.
         - garnetInCoeffs, coefficients controlling the P and T of the garnet in reaction
         - spinelOutCoeffs, coefficients controlling the P and T of the spinel out reaction
         - plagioclaseInInterval, The plagioclase in interval (in km).
        """
        # Check if using defaults, and assemble args if so:
        if method == 'default':
            default_kwargs = {'olv_D': pyMelt.chemistry.olv_D,
                              'cpx_D': pyMelt.chemistry.cpx_D,
                              'opx_D': pyMelt.chemistry.opx_D,
                              'spn_D': pyMelt.chemistry.spn_D,
                              'grt_D': pyMelt.chemistry.grt_D,
                              'plg_D': pyMelt.chemistry.plg_D,
                              'D': pyMelt.chemistry.workman05_D}

            for argname in default_kwargs:
                if argname not in kwargs:
                    kwargs[argname] = default_kwargs[argname]

            method = pyMelt.chemistry.default_methods

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
                elements = {self.mantle.names[0]: pyMelt.chemistry.workman05_dmm}
            elif elements is None:
                raise InputError("Either species_objects or elements must be provided.")

            species_objects = {}
            for lith in elements:
                method_recon = {}
                kwargs_recon = {}
                if(isinstance(method, dict)
                   and any(item in method.keys() for item in elements.keys())):
                    method_recon = method[lith]
                else:
                    method_recon = method
                for kw in kwargs:
                    if(isinstance(kwargs[kw], dict)
                       and any(item in kwargs[kw].keys() for item in elements.keys())):
                        kwargs_recon[kw] = kwargs[kw][lith]
                    else:
                        kwargs_recon[kw] = kwargs[kw]
                species_objects[lith] = self._create_species_objects(elements[lith],
                                                                     method_recon,
                                                                     **kwargs_recon)

        self._species_calc_type = {}
        for lith in species_objects:
            species_calc_type = []
            for species in species_objects[lith]:
                species_calc_type.append(species.calculation_type)
            self._species_calc_type[lith] = species_calc_type

        # Check that the lithology names are correct
        for lith in species_objects:
            if lith not in self.mantle.names:
                raise InputError("The lithology specified (" + lith + ") was not found in "
                                 "the mantle from which this melting column was constructed.")

        # Iterate through calculations for each lithology:
        for lith in species_objects:
            species_names = []
            for j in range(len(species_objects[lith])):
                species_names.append(species_objects[lith][j].name)
            results = np.zeros([np.shape(self.P)[0], len(species_objects[lith])])
            for i, row in self.lithologies[lith].iterrows():
                for j in range(len(species_objects[lith])):
                    if row.F > 1e-15:
                        results[i, j] = species_objects[lith][j].composition(row)
                    else:
                        results[i, j] = np.nan

            constructdf = pd.DataFrame(results, columns=species_names)

            # Checks if the element exists already, and drops in original if so:
            repeats = [value for value in species_names if value in self.lithologies[lith].columns]
            self.lithologies[lith].drop(repeats, inplace=True, axis=1)

            self.lithologies[lith] = pd.concat([self.lithologies[lith], constructdf], axis=1)

    def _create_species_objects(self, elements, method, **kwargs):
        methods = {'batch': pyMelt.chemistry.batchSpecies,
                   'continuous_instantaneous': pyMelt.chemistry.continuousSpecies_instantaneous,
                   'continuous_accumulated': pyMelt.chemistry.continuousSpecies_accumulated,
                   'invmel': pyMelt.chemistry.invmelSpecies}
        species_objects = []
        for el in elements:
            kwargs_recon = {}
            for arg in kwargs:
                if(isinstance(kwargs[arg], dict)
                   and any(item in kwargs[arg].keys() for item in elements.keys())):
                    if el in kwargs[arg].keys():
                        kwargs_recon[arg] = kwargs[arg][el]
                    else:
                        kwargs_recon[arg] = None
                else:
                    kwargs_recon[arg] = kwargs[arg]
            if isinstance(method, dict) and any(item in method.keys() for item in elements.keys()):
                species_objects.append(methods[method[el]](el, elements[el], **kwargs_recon))
            else:
                species_objects.append(methods[method](el, elements[el], **kwargs_recon))
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
