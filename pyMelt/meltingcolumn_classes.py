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


class MeltingColumn():
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

        self.lithologies = {}
        for i in range(self.mantle.number_lithologies):
            df = pd.DataFrame()
            df['Pressure'] = self.P
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

    def calculate_chemistry(self, elements=None, species_objects=None, method='batch',
                            output_type='accumulated', **kwargs):
        """
        Calculate the composition of the melt for each species object supplied. If species objects
        are not defined manually the default options can be used for each element by specifying
        them in the elements list. See documentation for more information about the default values.
        In built methods include batch melting, near-fractional melting, and INVMEL melting.

        Note that the effects of the aluminous phase transitions are not incorporated into the
        batch and near-fractional melting calculations.

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
        method : string
            One of 'batch', 'continuous', 'invmel'

        Notes
        -----

        The 'batch' melting routine uses:
         - D, the bulk partition coefficient

        The 'continuous' melting routine uses:
         - D, the bulk partition coefficient
         - phi, the porosity during melting (as a fraction, not percent).
        """

        # Assemble the species_objects dictionary if not provided
        if species_objects is None:
            if elements is None:
                raise InputError("Either species_objects or elements must be provided.")
            else:
                species_objects = {}
                for lith in elements:
                    kwargs_recon = {}
                    for kw in kwargs:
                        if isinstance(kwargs[kw], dict) and kwargs[kw].keys() == elements.keys():
                            kwargs_recon[kw] = kwargs[kw][lith]
                        else:
                            kwargs_recon[kw] = kwargs[kw]
                    species_objects[lith] = self._create_species_objects(elements[lith],
                                                                         method,
                                                                         **kwargs_recon)

        # Check that the lithology names are correct
        for lith in species_objects:
            if lith not in self.mantle.names:
                raise InputError("The lithology specified (" + lith + ") was not found in "
                                 "the mantle from which this melting column was constructed.")

        # Check that there isn't any mutually inconsistent chemistry happening.
        if self.chemistry_output is None:
            self.chemistry_output = output_type
        elif self.chemistry_output != output_type:
            raise InputError("The output_type must match the output_type used for previous "
                             "calculations with this melting column. To compare different "
                             "output formats create a duplicate column before calculating "
                             "chemistry.")

        # Iterate through calculations for each lithology:
        for lith in species_objects:
            species_names = []
            for j in range(len(species_objects[lith])):
                species_names.append(species_objects[lith][j].name)
            results = np.zeros([np.shape(self.P)[0], len(species_objects[lith])])
            for i, row in self.lithologies[lith].iterrows():
                for j in range(len(species_objects[lith])):
                    if row.F > 1e-15:
                        if output_type == 'accumulated':
                            results[i, j] = species_objects[lith][j].accumulated_melt(row)
                        elif output_type == 'instantaneous':
                            results[i, j] = species_objects[lith][j].instantaneous_melt(row)
                        else:
                            raise InputError("The output type was not recognised. It must be "
                                             "one of 'accumulated' or 'fractional'.")
                    else:
                        results[i, j] = np.nan

            constructdf = pd.DataFrame(results, columns=species_names)

            # Checks if the element exists already, and drops in original if so:
            repeats = [value for value in species_names if value in self.lithologies[lith].columns]
            self.lithologies[lith].drop(repeats, inplace=True, axis=1)

            self.lithologies[lith] = pd.concat([self.lithologies[lith], constructdf], axis=1)

    def _create_species_objects(self, elements, method, **kwargs):
        methods = {'batch': pyMelt.chemistry.BatchSpecies,
                   'continuous': pyMelt.chemistry.ContinuousSpecies,
                   'invmel': pyMelt.chemistry.invmelSpecies}
        species_objects = []
        for el in elements:
            kwargs_recon = {}
            for arg in kwargs:
                if isinstance(kwargs[arg], dict) and kwargs[arg].keys() == elements.keys():
                    kwargs_recon[arg] = kwargs[arg][el]
                else:
                    kwargs_recon[arg] = kwargs[arg]
            species_objects.append(methods[method](el, elements[el], **kwargs_recon))
        return species_objects
