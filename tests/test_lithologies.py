"""
Lithology Tests

This module provides testing routines for the pure lithology classes.
"""

import unittest
import pyMelt as m
import numpy as np



class test_lithology(unittest.TestCase):
    def setUp(self):
        self.lithologies_to_test = {
                                     'katz_lherzolite': m.lithologies.katz.lherzolite(),
                                   # 'matthews_klb1': m.lithologies.matthews.klb1(),
                                   # 'matthews_kg1': m.lithologies.matthews.kg1(),
                                   # 'matthews_eclogite': m.lithologies.matthews.eclogite(),
                                   # 'pertermann_g2': m.lithologies.pertermann.g2(),
                                   # 'shorttle_kg1': m.lithologies.shorttle.kg1(),
                                   # 'shorttle_harzburgite': m.lithologies.shorttle.harzburgite(),
                                     }

        self.Tsolidus_ref = {
                             'katz_lherzolite': np.array([1085.7, 1213.5, 1331.1, 1438.5, 1535.7,
                                                          1622.7, 1699.5, 1766.1, 1822.5]),
                             }

        self.Tliquidus_ref = {
                              'katz_lherzolite': np.array([1780., 1823., 1862., 1897., 1928.,
                                                           1955., 1978., 1997., 2012.]),
                              }

        self.F_ref = {'katz_lherzolite': np.array([0.        , 0.        , 0.        , 0.        ,
                                                   0.        , 0.        , 0.        , 0.12371537,
                                                   0.27606232, 0.46115328, 0.71020361, 1.        ,
                                                   1.        ])
                      }

        self.dTdF_ref = {'katz_lherzolite': np.array([     np.inf,       np.inf,       np.inf,
                                                           np.inf,       np.inf,       np.inf,
                                                           np.inf, 331.40586724, 381.92595966,
                                                     321.8842976 , 278.73227917,      np.inf,
                                                           np.inf])
                         }

        self.dTdP_ref = {'katz_lherzolite': np.array([1.21212121e-02, 1.21212121e-02,
                                                      1.21212121e-02, 1.21212121e-02,
                                                      1.21212121e-02, 1.21212121e-02,
                                                      1.21212121e-02, 9.19962051e+01,
                                                      8.56167104e+01, 6.55068107e+01,
                                                      4.85828592e+01, 2.34482759e-02,
                                                      2.34482759e-02])
                         }

        self.P = np.linspace(0, 8, 9)
        self.T = np.linspace(800.0, 2000.0, 13)

    def runTest(self):
        for lith_name in list(self.lithologies_to_test.keys()):
            T = self.lithologies_to_test[lith_name].TSolidus(self.P)
            for i in range(len(self.P)):
                self.assertAlmostEqual(T[i], self.Tsolidus_ref[lith_name][i],
                                       msg=lith_name + ' solidus')

            T = self.lithologies_to_test[lith_name].TLiquidus(self.P)
            for i in range(len(self.P)):
                self.assertAlmostEqual(T[i], self.Tliquidus_ref[lith_name][i],
                                       msg=lith_name + ' liquidus')

            for i in range(len(self.T)):
                F = self.lithologies_to_test[lith_name].F(P=3.0, T=self.T[i])
                self.assertAlmostEqual(F, self.F_ref[lith_name][i], places=7, msg=lith_name+' F')

            for i in range(len(self.T)):
                dTdF = self.lithologies_to_test[lith_name].dTdF(P=3.0, T=self.T[i])
                self.assertAlmostEqual(dTdF, self.dTdF_ref[lith_name][i], places=7,
                                       msg=lith_name+' dTdF')

            for i in range(len(self.T)):
                dTdP = self.lithologies_to_test[lith_name].dTdP(P=3.0, T=self.T[i])
                self.assertAlmostEqual(dTdP, self.dTdP_ref[lith_name][i], places=7,
                                       msg=lith_name+' dTdP')
