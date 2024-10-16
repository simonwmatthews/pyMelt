"""
Chemistry Tests

This module provides testing routines for the Chemistry classes.
"""

import pytest
import unittest
import pyMelt as m
from numpy import allclose
import numpy as np

def test_should_calculate_majors_with_isobaric_melting_start():
    lz = m.lithologies.matthews.klb1()
    mantle = m.mantle([lz], [1.0], ['lz'])
    column = mantle.adiabaticMelt(2000.0, Pstart=8.0)
    column.calculateMajorOxides()
    target = 36.338872220258374
    testcomp = column.composition['lz'].liq_MgO.iloc[0]
    assert allclose(target, testcomp)

def test_should_calculate_traces_with_isobaric_melting_start():
    lz = m.lithologies.matthews.klb1()
    mantle = m.mantle([lz], [1.0], ['lz'])
    column = mantle.adiabaticMelt(1400.0, Pstart=8.0)
    column.calculateMineralProportions()
    column.calculateTraceElements(c0={'lz':{'La':1.0}})
    F = np.array(column.composition['lz'].F)
    dF = F[1:] - F[:-1]
    c = np.array(column.composition['lz'].liq_La)[1:]
    total = np.sum(c[~np.isnan(c)] * dF[~np.isnan(c)])
    print(total)
    print(column.composition['lz'].liq_La)
    print(column.F)
    assert False

# These are old tests for the previous version of pyMelt and need to be rewritten for the new chemistry interface
@pytest.mark.skip
class test_chemistry(unittest.TestCase):

    def setUp(self):
        self.lz = m.lithologies.matthews.klb1()
        self.px = m.lithologies.matthews.kg1()

    def test_single_lithology_column(self):
        """
        Test column chemistry is correctly calculated for La, Dy, and Yb during
        single lithology melting.
        """
        mantle = m.mantle([self.lz], [1.0], ['lz'])
        column = mantle.adiabaticMelt(1350.0)
        column.calculateChemistry()

        concentrations = {'La': 26.029913793159484,
                          'Dy': 9.613940884597207,
                          'Yb': 4.497467642416129}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(column.lithologies['lz'][element].iloc[1],
                                   concentrations[element], places=6,
                                   msg="Single lithology column " + element + " is incorrect.")

    def test_single_lithology_spreadingCentre(self):
        """
        Test column chemistry is correctly calculated for La, Dy, and Yb during
        single lithology melting and spreading centre homogenisation.
        """
        mantle = m.mantle([self.lz], [1.0], ['lz'])
        column = mantle.adiabaticMelt(1350.0)
        column.calculateChemistry()
        morb = m.geosettings.spreadingCentre(column)

        concentrations = {'La': 1.5115734212608263,
                          'Dy': 3.37049583704782,
                          'Yb': 2.1831968869748395}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(morb.chemistry[element],
                                   concentrations[element], places=6,
                                   msg=("Single lithology homogenised sc " + element
                                        + " is incorrect."))

    def test_single_lithology_intraPlate(self):
        """
        Test column chemistry is correctly calculated for La, Dy, and Yb during
        single lithology melting and intraplate homogenisation.
        """
        mantle = m.mantle([self.lz], [1.0], ['lz'])
        column = mantle.adiabaticMelt(1350.0)
        column.calculateChemistry()
        oib = m.geosettings.intraPlate(column, 1.0,
                                       weightingFunction=m.geosettings.weighting_expdecay,
                                       weighting_wavelength=0.1)

        concentrations = {'La': 7.247433001514023,
                          'Dy': 7.339268956958773,
                          'Yb': 3.764364568503429}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(oib.chemistry[element],
                                   concentrations[element], places=6,
                                   msg=("Single lithology homogenised ip " + element
                                        + " is incorrect."))

    def test_multi_lithology_column(self):
        """
        Test column chemistry is correctly calculated for La, Dy, and Yb during
        multi-lithology melting..
        """
        mantle = m.mantle([self.lz, self.px], [1.0, 1.0], ['lz', 'px'])
        column = mantle.adiabaticMelt(1350.0)
        column.calculateChemistry(elements={'lz': m.chemistry.workman05_dmm,
                                            'px': m.chemistry.stracke03_bsic},
                                  cpxExhaustion={'lz': 0.18,
                                                 'px': 0.70},
                                  garnetInCoeffs={'lz': [666.7, 400.0],
                                                  'px': [666.7, 400.0]},
                                  spinelOutCoeffs={'lz': [666.7, 533.0],
                                                   'px': [666.7, 533.0]},
                                  mineralProportions={'lz': m.chemistry.klb1_MineralProportions,
                                                      'px': m.chemistry.kg1_MineralProportions}
                                  )

        lz_conc = {'La': 0.41640981547719697,
                   'Dy': 2.486516755289274,
                   'Yb': 2.271169662741649}

        px_conc = {'La': 78.64457868423119,
                   'Dy': 4.92944558747092,
                   'Yb': 1.101745291927927}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(column.lithologies['lz'][element].iloc[-1],
                                   lz_conc[element], places=6,
                                   msg="Multi-lithology column lz " + element + " is incorrect.")

            self.assertAlmostEqual(column.lithologies['px'][element].iloc[1],
                                   px_conc[element], places=6,
                                   msg="Multi-lithology column px " + element + " is incorrect.")

    def test_multi_lithology_spreadingCentre(self):
        """
        Test column chemistry is correctly calculated for La, Dy, and Yb during
        multi-lithology melting and spreading centre homogenisation.
        """
        mantle = m.mantle([self.lz, self.px], [1.0, 1.0], ['lz', 'px'])
        column = mantle.adiabaticMelt(1350.0)
        column.calculateChemistry(elements={'lz': m.chemistry.workman05_dmm,
                                            'px': m.chemistry.stracke03_bsic},
                                  cpxExhaustion={'lz': 0.18,
                                                 'px': 0.70},
                                  garnetInCoeffs={'lz': [666.7, 400.0],
                                                  'px': [666.7, 400.0]},
                                  spinelOutCoeffs={'lz': [666.7, 533.0],
                                                   'px': [666.7, 533.0]},
                                  mineralProportions={'lz': m.chemistry.klb1_MineralProportions,
                                                      'px': m.chemistry.kg1_MineralProportions}
                                  )

        morb = m.geosettings.spreadingCentre(column)

        concentrations = {'La': 4.281535100966795,
                          'Dy': 9.692041695712126,
                          'Yb': 5.603085296135785}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(morb.chemistry[element],
                                   concentrations[element], places=6,
                                   msg=("Multi-lithology homogenised sc " + element
                                        + " is incorrect."))

    def test_multi_lithology_intraPlate(self):
        """
        Test column chemistry is correctly calculated for La, Dy, and Yb during
        multi-lithology melting and intraplate homogenisation.
        """
        mantle = m.mantle([self.lz, self.px], [1.0, 1.0], ['lz', 'px'])
        column = mantle.adiabaticMelt(1350.0)
        column.calculateChemistry(elements={'lz': m.chemistry.workman05_dmm,
                                            'px': m.chemistry.stracke03_bsic},
                                  cpxExhaustion={'lz': 0.18,
                                                 'px': 0.70},
                                  garnetInCoeffs={'lz': [666.7, 400.0],
                                                  'px': [666.7, 400.0]},
                                  spinelOutCoeffs={'lz': [666.7, 533.0],
                                                   'px': [666.7, 533.0]},
                                  mineralProportions={'lz': m.chemistry.klb1_MineralProportions,
                                                      'px': m.chemistry.kg1_MineralProportions}
                                  )

        oib = m.geosettings.intraPlate(column, 1.0,
                                       weightingFunction=m.geosettings.weighting_expdecay,
                                       weighting_wavelength=0.1)

        concentrations = {'La': 22.731167287563885,
                          'Dy': 10.841739275444112,
                          'Yb': 5.126314925222699}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(oib.chemistry[element],
                                   concentrations[element], places=6,
                                   msg=("Multi-lithology homogenised ip " + element
                                        + " is incorrect."))

    def test_variable_D(self):
        """
        Test that distribution coefficients can be changed for all cases.
        """
        mantle = m.mantle([self.lz, self.px], [1.0, 1.0], ['lz', 'px'])
        column = mantle.adiabaticMelt(1350.0)

        olv_D_new = m.chemistry.olv_D
        olv_D_new['La'] = 1.3
        olv_D_new['Dy'] = 1.3
        olv_D_new['Yb'] = 1.3

        column.calculateChemistry(olv_D=olv_D_new,
                                  elements={'lz': m.chemistry.workman05_dmm,
                                            'px': m.chemistry.stracke03_bsic},
                                  cpxExhaustion={'lz': 0.18,
                                                 'px': 0.70},
                                  garnetInCoeffs={'lz': [666.7, 400.0],
                                                  'px': [666.7, 400.0]},
                                  spinelOutCoeffs={'lz': [666.7, 533.0],
                                                   'px': [666.7, 533.0]},
                                  mineralProportions={'lz': m.chemistry.klb1_MineralProportions,
                                                      'px': m.chemistry.kg1_MineralProportions}
                                  )

        lz_conc_D_new = {'La': 0.2110498664868391,
                         'Dy': 0.5497607395120011,
                         'Yb': 0.39177925047563783}

        px_conc_D_new = {'La': 6.547065097730284,
                         'Dy': 4.004202514670137,
                         'Yb': 1.0150774093931085}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(column.lithologies['lz'][element].iloc[-1],
                                   lz_conc_D_new[element], places=6,
                                   msg=("Multi-lithology column lz " + element
                                        + " (D) is incorrect."))

            self.assertAlmostEqual(column.lithologies['px'][element].iloc[1],
                                   px_conc_D_new[element], places=6,
                                   msg=("Multi-lithology column px " + element
                                        + " (D) is incorrect."))

        morb_new = m.geosettings.spreadingCentre(column)

        concentrations_spreadingCentre = {'La': 3.3944776138412607,
                                          'Dy': 7.8304433811896805,
                                          'Yb': 4.520023486233671}

        oib_new = m.geosettings.intraPlate(column, 1.0,
                                           weightingFunction=m.geosettings.weighting_expdecay,
                                           weighting_wavelength=0.1)

        concentrations_intraPlate = {'La': 6.296271631975634,
                                     'Dy': 7.725888101870476,
                                     'Yb': 3.7066095658285123}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(morb_new.chemistry[element],
                                   concentrations_spreadingCentre[element], places=6,
                                   msg=("Multi-lithology homogenised sc " + element
                                        + " (D) is incorrect."))

            self.assertAlmostEqual(oib_new.chemistry[element],
                                   concentrations_intraPlate[element], places=6,
                                   msg=("Multi-lithology homogenised ip " + element
                                        + " (D) is incorrect."))
