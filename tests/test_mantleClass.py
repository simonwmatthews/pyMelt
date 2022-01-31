"""
Mantle Class Tests

This module provides testing routines for the Mantle class
"""

import unittest
import pyMelt as m
import numpy as np

from pyMelt.core import InputError


class test_mantle(unittest.TestCase):
    def setUp(self):
        self.models = {'katz_lherzolite': m.lithologies.katz.lherzolite(),
                       'matthews_klb1': m.lithologies.matthews.klb1(),
                       'matthews_kg1': m.lithologies.matthews.kg1(),
                       'matthews_eclogite': m.lithologies.matthews.eclogite(),
                       'mckenzie_lherzolite': m.lithologies.mckenzie.lherzolite(),
                       'pertermann_g2': m.lithologies.pertermann.g2(),
                       'shorttle_kg1': m.lithologies.shorttle.kg1()
                       }

    def test_single_lithology(self):
        """
        Test adiabtic decompression in pure lithologies using default calculation settings.
        """

        F_top_1300 = {'katz_lherzolite': 0.1840275294566579,
                      'matthews_klb1': 0.19171243561940424,
                      'matthews_kg1': 0.4197494916507786,
                      'matthews_eclogite': 0.39128784694385865,
                      'mckenzie_lherzolite': 0.25189890081706984,
                      'pertermann_g2': 0.46778134851714903,
                      'shorttle_kg1': 0.30985171334907324
                      }

        T_top_1300 = {'katz_lherzolite': 1213.1669678349292,
                      'matthews_klb1': 1210.5333533430478,
                      'matthews_kg1': 1107.8615259442008,
                      'matthews_eclogite': 1116.0376787493344,
                      'mckenzie_lherzolite': 1211.5159620032775,
                      'pertermann_g2': 1083.8813041535018,
                      'shorttle_kg1': 1141.927174068584
                      }

        F_top_1500 = {'katz_lherzolite': 0.3482562827506706,
                      'matthews_klb1': 0.33831997369443717,
                      'matthews_kg1': 0.653972152744598,
                      'matthews_eclogite': 0.6298137735462399,
                      'mckenzie_lherzolite': 0.36672906306884334,
                      'pertermann_g2': 0.7146561815636401,
                      'shorttle_kg1': 0.5838763866298998
                      }

        T_top_1500 = {'katz_lherzolite': 1305.3033208407003,
                      'matthews_klb1': 1331.6721589237834,
                      'matthews_kg1': 1163.088871244221,
                      'matthews_eclogite': 1161.9929806182686,
                      'mckenzie_lherzolite': 1348.8539471940846,
                      'pertermann_g2': 1131.185304594791,
                      'shorttle_kg1': 1168.7766254424896
                      }

        for lith in self.models.keys():
            mantle = m.mantle([self.models[lith]], [1.0])
            col1300 = mantle.adiabaticMelt(1300.0)
            col1500 = mantle.adiabaticMelt(1500.0)

            self.assertAlmostEqual(col1300.F.iloc[-1], F_top_1300[lith], places=6,
                                   msg=lith + 'F at Tp=1300')
            self.assertAlmostEqual(col1500.F.iloc[-1], F_top_1500[lith], places=6,
                                   msg=lith + 'F at Tp=1500')

            self.assertAlmostEqual(col1300.T.iloc[-1], T_top_1300[lith], places=6,
                                   msg=lith + 'T at Tp=1300')
            self.assertAlmostEqual(col1500.T.iloc[-1], T_top_1500[lith], places=6,
                                   msg=lith + 'T at Tp=1500')

    def test_melt_steps(self):
        """
        Test that steps can be calculated
        """
        mantle = m.mantle([self.models['matthews_klb1']], [1.0])

        col = mantle.adiabaticMelt(1300, steps=100)
        self.assertEqual(np.shape(col.F)[0], 100,
                         msg="Incorrect number of steps")

        self.assertAlmostEqual(col.P.iloc[0], 1.403882594274345, places=6,
                               msg="Incorrect starting pressure.")

        self.assertAlmostEqual(col.P.iloc[-1], 0.01, places=6,
                               msg="Incorrect final pressure.")

    def test_melt_steps_startP(self):
        """
        Test that steps with a starting pressure can be calculated.
        """
        mantle = m.mantle([self.models['matthews_klb1']], [1.0])

        col = mantle.adiabaticMelt(1300, steps=100, Pstart=5.0)
        self.assertEqual(np.shape(col.F)[0], 100,
                         msg="Incorrect number of steps when Pstart specified.")

        self.assertAlmostEqual(col.P.iloc[0], 4.9825594629612135, places=6,
                               msg="Incorrectly adjusted starting pressure.")

        self.assertEqual(col.P.iloc[-1], 0.0, msg="Incorrectly adjusted final pressure.")

    def test_melt_steps_startP_noAdj(self):
        """
        Test that steps with a starting pressure can be calculated without adjustment
        """
        mantle = m.mantle([self.models['matthews_klb1']], [1.0])

        col = mantle.adiabaticMelt(1300, steps=100, Pstart=5.0, adjust_pressure=False, Pend=1.0)
        self.assertEqual(np.shape(col.F)[0], 100,
                         msg="Incorrect number of steps when Pstart specified.")

        self.assertEqual(col.P.iloc[0], 5.0, msg="Starting pressure not propagated.")

        self.assertEqual(col.P.iloc[-1], 1.0, msg="Final pressure not propagated.")

    def test_nonmelting_auto(self):
        """
        Test that a non-melting mantle will return an InputError
        """
        hz = m.lithologies.shorttle.harzburgite()
        mantle = m.mantle([hz], [1.0])
        self.assertRaises(InputError, mantle.adiabaticMelt, 1300.0)

    def test_nonmelting_Prange(self):
        """
        Test that a non-melting mantle can produce an adiabat when pressure is specified.
        """
        hz = m.lithologies.shorttle.harzburgite()
        mantle = m.mantle([hz], [1.0])
        col = mantle.adiabaticMelt(1300, Pstart=5.0)

        self.assertAlmostEqual(col.T.iloc[-1], 1300, places=0,
                               msg="Tp not returned by non-melting adiabat.")

    def test_adiabat(self):
        """
        Test that the adiabat is calculated correctly. Tests two lithologies only, just to check
        the different constants are propagated properly.
        """
        hz = m.lithologies.shorttle.harzburgite()
        mantle = m.mantle([hz], [1.0])
        self.assertAlmostEqual(mantle.adiabat(5.0, 1300.0), 1374.308545262, places=6,
                               msg="Adiabat temperature for Shorttle-Harzburgite is not "
                                   "predicted correctly.")

        mantle = m.mantle([self.models['matthews_klb1']], [1.0])
        self.assertAlmostEqual(mantle.adiabat(5.0, 1300), 1398.2908507468187, places=6,
                               msg="Adiabat temperature for Matthews-KLB1 is not "
                                   "predicted correctly.")

    def test_superSolidusStart(self):
        """
        Check that a super solidus start is identified and calculated correctly.
        """
        mantle = m.mantle([self.models['matthews_eclogite']], [1.0])

        self.assertWarns(UserWarning, mantle.adiabaticMelt, Tp=1500, Pstart=3.0)

        col = mantle.adiabaticMelt(Tp=1500, Pstart=3.0, ReportSSS=False)

        self.assertAlmostEqual(col.T.iloc[0], 1393.4183401385737, places=6,
                               msg="Isobaric melting stop not working correctly.")

    def test_multilithology_auto(self):
        """
        Test that multi-lithology melting is calculated correctly.
        """
        mantle = m.mantle([self.models['matthews_klb1'], self.models['matthews_kg1']],
                          [1.0, 0.25],
                          ['klb1', 'kg1'])

        col = mantle.adiabaticMelt(1350.0)

        self.assertAlmostEqual(col.F.iloc[-1], 0.3031835235468895, places=6,
                               msg="Total melt fraction is incorrect.")

        self.assertAlmostEqual(col.lithologies['klb1'].F.iloc[-1], 0.16857773509268661, places=6,
                               msg="KLB1 melt fraction is incorrect.")

        self.assertAlmostEqual(col.lithologies['kg1'].F.iloc[-1], 0.8416066773637008, places=6,
                               msg="KG1 melt fraction is incorrect.")

        self.assertAlmostEqual(col.T.iloc[-1], 1204.0000266711634, places=6,
                               msg="The temperature at the top of the column is incorrect.")

    def test_prevent_freezing(self):
        """
        Test that pyMelt will prevent freezing
        """
        mantle = m.mantle([self.models['matthews_klb1'],
                           self.models['matthews_kg1']],
                          [10.0, 1.0],
                          ['klb1', 'kg1'])

        self.assertWarns(UserWarning, mantle.adiabaticMelt, Tp=1600.0)

        col = mantle.adiabaticMelt(1600.0, warn_prevent_freezing=False)
        col2 = mantle.adiabaticMelt(1600.0, prevent_freezing=False)

        s1 = np.shape(col.lithologies['klb1'].F[col.lithologies['klb1'].F > 0].unique())[0]
        s2 = np.shape(col2.lithologies['klb1'].F[col2.lithologies['klb1'].F > 0].unique())[0]

        self.assertTrue(s1 < s2, msg="No evidence of freezing prevention")
