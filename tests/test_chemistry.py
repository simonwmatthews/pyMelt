"""
Chemistry Tests

This module provides testing routines for the Chemistry classes.
"""

import unittest
import pyMelt as m


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

        concentrations = {'La': 25.966120238876435,
                          'Dy': 9.611454235124366,
                          'Yb': 4.496720428496521}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(column.lithologies['lz'][element].iloc[1],
                                   concentrations[element], places=6,
                                   msg="Single lithology column "+element+" is incorrect.")

    def test_single_lithology_spreadingCentre(self):
        """
        Test column chemistry is correctly calculated for La, Dy, and Yb during
        single lithology melting and spreading centre homogenisation.
        """
        mantle = m.mantle([self.lz], [1.0], ['lz'])
        column = mantle.adiabaticMelt(1350.0)
        column.calculateChemistry()
        morb = m.geosettings.spreadingCentre(column)

        concentrations = {'La': 1.4722236338647683,
                          'Dy': 2.8617962509280677,
                          'Yb': 1.97280007634921}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(morb.chemistry[element],
                                   concentrations[element], places=6,
                                   msg="Single lithology homogenised sc "+element+" is incorrect.")

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

        concentrations = {'La': 7.104201514862648,
                          'Dy': 6.879161621049682,
                          'Yb': 3.611603816444381}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(oib.chemistry[element],
                                   concentrations[element], places=6,
                                   msg="Single lithology homogenised ip "+element+" is incorrect.")

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

        lz_conc = {'La': 0.29052786488937843,
                   'Dy': 1.7951039754820055,
                   'Yb': 1.9519669624938756}

        px_conc = {'La': 78.5201430023945,
                   'Dy': 4.929521209273142,
                   'Yb': 1.1017855590754706}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(column.lithologies['lz'][element].iloc[-1],
                                   lz_conc[element], places=6,
                                   msg="Multi-lithology column lz "+element+" is incorrect.")

            self.assertAlmostEqual(column.lithologies['px'][element].iloc[1],
                                   px_conc[element], places=6,
                                   msg="Multi-lithology column px "+element+" is incorrect.")

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

        concentrations = {'La': 4.235660539351117,
                          'Dy': 8.703630112867264,
                          'Yb': 5.091522033090427}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(morb.chemistry[element],
                                   concentrations[element], places=6,
                                   msg="Multi-lithology homogenised sc "+element+" is incorrect.")

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

        concentrations = {'La': 22.57146561994308,
                          'Dy': 10.474179181111653,
                          'Yb': 4.954627191741076}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(oib.chemistry[element],
                                   concentrations[element], places=6,
                                   msg="Multi-lithology homogenised ip "+element+" is incorrect.")

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

        lz_conc_D_new = {'La': 0.22039527545369597,
                         'Dy': 0.571489745287995,
                         'Yb': 0.4076187542115444}

        px_conc_D_new = {'La': 6.546043921246597,
                         'Dy': 4.004231815822892,
                         'Yb': 1.0151094027031582}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(column.lithologies['lz'][element].iloc[-1],
                                   lz_conc_D_new[element], places=6,
                                   msg="Multi-lithology column lz "+element+" (D) is incorrect.")

            self.assertAlmostEqual(column.lithologies['px'][element].iloc[1],
                                   px_conc_D_new[element], places=6,
                                   msg="Multi-lithology column px "+element+" (D) is incorrect.")

        morb_new = m.geosettings.spreadingCentre(column)

        concentrations_spreadingCentre = {'La': 3.878965128767849,
                                          'Dy': 8.118681406036563,
                                          'Yb': 4.735525468396869}

        oib_new = m.geosettings.intraPlate(column, 1.0,
                                           weightingFunction=m.geosettings.weighting_expdecay,
                                           weighting_wavelength=0.1)

        concentrations_intraPlate = {'La': 6.5793025139526184,
                                     'Dy': 7.777604213908181,
                                     'Yb': 3.7422394633859084}

        for element in ['La', 'Dy', 'Yb']:
            self.assertAlmostEqual(morb_new.chemistry[element],
                                   concentrations_spreadingCentre[element], places=6,
                                   msg="Multi-lithology homogenised sc "+element+" (D) is incorrect.")

            self.assertAlmostEqual(oib_new.chemistry[element],
                                   concentrations_intraPlate[element], places=6,
                                   msg="Multi-lithology homogenised ip "+element+" (D) is incorrect.")
