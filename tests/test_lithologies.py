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
                                    'matthews_klb1': m.lithologies.matthews.klb1(),
                                    'matthews_kg1': m.lithologies.matthews.kg1(),
                                    'matthews_eclogite': m.lithologies.matthews.eclogite(),
                                   'pertermann_g2': m.lithologies.pertermann.g2(),
                                   'shorttle_kg1': m.lithologies.shorttle.kg1(),
                                   # 'shorttle_harzburgite': m.lithologies.shorttle.harzburgite(),
                                     }

        self.Tsolidus_ref = {
                             'katz_lherzolite': np.array([1085.7, 1213.5, 1331.1, 1438.5, 1535.7,
                                                          1622.7, 1699.5, 1766.1, 1822.5]),
                             'matthews_klb1': np.array([1130.35520902, 1275.08331801,
                                                       1397.57332987, 1501.53525181,
                                                       1589.82173407, 1664.6734521 ,
                                                       1727.88253266, 1780.90482445,
                                                       1824.93910146]),
                             'matthews_kg1': np.array([ 957.27102949, 1149.6685339 , 1292.55257242,
                                                       1407.80973768, 1505.41038185, 1590.7438744 ,
                                                       1667.05570292, 1736.45257412, 1800.379823  ]
                                                       ),
                             'matthews_eclogite': np.array([ 931.56250599, 1050.46785495,
                                                            1153.92452251, 1246.11834014,
                                                            1329.73512912, 1406.60153262,
                                                            1478.01641838, 1544.93645645,
                                                            1608.08668088]),
                             'pertermann_g2': np.array([ 920., 1050., 1180., 1310., 1440., 1570.,
                                                        1700., 1830., 1960.]),
                             'shorttle_kg1': np.array([1095.4, 1214.8, 1324.8, 1425.4, 1516.6,
                                                       1598.4, 1670.8, 1733.8, 1787.4]),
                             'mckenzie_lz': np.array([1099.9324899135006, 1235.6558912187331,
                                                      1370.2692364054583, 1499.8089711536893,
                                                      1612.3776612480988, 1694.8871653527426,
                                                      1750.3691421541535, 1788.9254622582778,
                                                      1817.4674524991735])
                             }

        self.Tliquidus_ref = {
                              'katz_lherzolite': np.array([1780., 1823., 1862., 1897., 1928.,
                                                           1955., 1978., 1997., 2012.]),
                              'matthews_klb1': np.array([1885.29107647, 1898.28001541,
                                                         1911.26789493, 1924.25471815,
                                                         1937.24048821, 1950.22520822,
                                                         1963.20888128, 1976.19151047,
                                                         1989.17309886]),
                              'matthews_kg1': np.array([1519.00306938, 1586.28255841,
                                                        1653.56051792, 1720.83695691,
                                                        1788.11188432, 1855.38530902,
                                                        1922.65723977, 1989.92768528,
                                                        2057.19665418]),
                              'matthews_eclogite': np.array([1216.0012043 , 1320.0604384 ,
                                                             1412.02002231, 1494.17031038,
                                                             1568.20720152, 1635.4215767 ,
                                                             1696.8188124 , 1753.19722384,
                                                             1805.20129058]),
                              'pertermann_g2': np.array([1175., 1289., 1403., 1517., 1631., 1745.,
                                                        1859., 1973., 2087.]),
                              'shorttle_kg1': np.array([1780., 1823., 1862., 1897., 1928., 1955.,
                                                        1978., 1997., 2012.]),
                              'mckenzie_lz': np.array([1736.        , 1816.61835502, 1876.80736135,
                                                       1917.21992504, 1945.02072556, 1965.33801856,
                                                       1981.09476852, 1993.93714821, 2004.82809674]
                                                       )
                              }

        self.F_ref = {'katz_lherzolite': np.array([0.        , 0.        , 0.        , 0.        ,
                                                   0.        , 0.        , 0.        , 0.12371537,
                                                   0.27606232, 0.46115328, 0.71020361, 1.        ,
                                                   1.        ]),
                      'matthews_klb1': np.array([0.        , 0.        , 0.        , 0.        ,
                                                 0.        , 0.        , 0.        , 0.        ,
                                                 0.24453539, 0.40844828, 0.64044998, 0.92432968,
                                                 1.        ]),
                      'matthews_kg1': np.array([0.        , 0.        , 0.        , 0.        ,
                                                0.        , 0.        , 0.        , 0.37883722,
                                                0.5474943 , 0.90446444,1.        , 1.        ,
                                                1.        ]),
                      'matthews_eclogite': np.array([0.        , 0.        , 0.        ,
                                                     0.        , 0.        , 0.03845395,
                                                     0.36099624, 1.        , 1.        ,
                                                     1.        , 1.        , 1.        ,
                                                     1.        ]),
                      'pertermann_g2': np.array([0.        , 0.        , 0.        , 0.        ,
                                                 0.        , 0.        , 0.25371645, 0.86233368,
                                                 1.        , 1.        , 1.        , 1.        ,
                                                 1.        ]),
                      'shorttle_kg1': np.array([0.        , 0.        , 0.        , 0.        ,
                                                0.        , 0.        , 0.        , 0.35806323,
                                                0.7482411 , 0.80917489, 0.89638302, 1.        ,
                                                1.        ])
                      }

        self.dTdF_ref = {'katz_lherzolite': np.array([     np.inf,       np.inf,       np.inf,
                                                           np.inf,       np.inf,       np.inf,
                                                           np.inf, 331.40586724, 381.92595966,
                                                     321.8842976 , 278.73227917,      np.inf,
                                                           np.inf]),
                        'matthews_klb1': np.array([      np.inf,       np.inf,       np.inf,
                                                         np.inf,       np.inf,       np.inf,
                                                         np.inf,       np.inf,  402.0305668,
                                                   338.83944065, 291.66069592, 258.08574308,
                                                         np.inf]),
                        'matthews_kg1': np.array([      np.inf,       np.inf,       np.inf,
                                                        np.inf,       np.inf,       np.inf,
                                                        np.inf, 202.82636818, 172.20500314,
                                                  137.76902666,       np.inf,       np.inf,
                                                        np.inf]),
                        'matthews_eclogite': np.array([      np.inf,       np.inf,       np.inf,
                                                             np.inf,       np.inf, 656.60715452,
                                                       199.75136891,       np.inf,       np.inf,
                                                             np.inf,       np.inf,       np.inf,
                                                             np.inf]),
                         'pertermann_g2': np.array([      np.inf,       np.inf,       np.inf,
                                                          np.inf,       np.inf,       np.inf,
                                                    229.00873514, 128.11152331,       np.inf,
                                                          np.inf,       np.inf,       np.inf,
                                                          np.inf]),
                         'shorttle_kg1': np.array([       np.inf,        np.inf,        np.inf,
                                                          np.inf,        np.inf,        np.inf,
                                                          np.inf,  158.72747585, 2310.02208704,
                                                   1321.98074238, 1022.21610823,        np.inf,
                                                          np.inf])

                         }

        self.dTdP_ref = {'katz_lherzolite': np.array([1.21212121e-02, 1.21212121e-02,
                                                      1.21212121e-02, 1.21212121e-02,
                                                      1.21212121e-02, 1.21212121e-02,
                                                      1.21212121e-02, 9.19962051e+01,
                                                      8.56167104e+01, 6.55068107e+01,
                                                      4.85828592e+01, 2.34482759e-02,
                                                      2.34482759e-02]),
                         'matthews_klb1': np.array([1.21212121e-02, 1.21212121e-02, 1.21212121e-02,
                                                    1.21212121e-02, 1.21212121e-02, 1.21212121e-02,
                                                    1.21212121e-02, 1.21212121e-02, 3.47790758e+01,
                                                    2.79486747e+01, 2.12218674e+01, 1.45841906e+01,
                                                    2.34482759e-02]),
                         'matthews_kg1': np.array([1.21212121e-02, 1.21212121e-02, 1.21212121e-02,
                                                   1.21212121e-02, 1.21212121e-02, 1.21212121e-02,
                                                   1.21212121e-02, 9.88283751e+01, 8.40880560e+01,
                                                   7.01057934e+01, 2.34482759e-02, 2.34482759e-02,
                                                   2.34482759e-02]),
                         'matthews_eclogite': np.array([1.21212121e-02, 1.21212121e-02,
                                                        1.21212121e-02, 1.21212121e-02,
                                                        1.21212121e-02, 8.54380136e+01,
                                                        8.15298753e+01, 2.34482759e-02,
                                                        2.34482759e-02, 2.34482759e-02,
                                                        2.34482759e-02, 2.34482759e-02,
                                                        2.34482759e-02]),
                         'pertermann_g2': np.array([1.21212121e-02, 1.21212121e-02, 1.21212121e-02,
                                                    1.21212121e-02, 1.21212121e-02, 1.21212121e-02,
                                                    1.23043478e+02, 1.15314010e+02, 1.21212121e-02,
                                                    1.21212121e-02, 1.21212121e-02, 1.21212121e-02,
                                                    1.21212121e-02]),
                         'shorttle_kg1': np.array([1.22227403e+02, 1.18017712e+02, 1.13808022e+02,
                                                   1.09598332e+02, 1.05388642e+02, 1.01178952e+02,
                                                   9.69692613e+01, 9.27595711e+01, 8.24856812e+01,
                                                   6.58238357e+01, 4.91619902e+01, 1.21212121e-02,
                                                   1.21212121e-02])
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
                self.assertAlmostEqual(dTdP, self.dTdP_ref[lith_name][i], places=6,
                                       msg=lith_name+' dTdP')
