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

    def test_single_lithology_auto(self):
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

    def test_nonmelting_auto(self):
        """
        Test that a non-melting mantle will return an InputError
        """
        hz = m.lithologies.shorttle.harzburgite()
        mantle = m.mantle([hz], [1.0])
        self.assertRaises(InputError, mantle.adiabaticMelt, 1300.0)
