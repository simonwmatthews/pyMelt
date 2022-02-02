"""
Hydrous Lithology Tests

This module provides testing routines for the pure lithology classes.
"""

import unittest
import pyMelt as m
import numpy as np


class test_lithologyProperties(unittest.TestCase):
    def setUp(self):
        lz = m.lithologies.matthews.klb1()
        self.hlz = m.hydrousLithology(lz, 1.0)
        self.P = np.linspace(0, 8, 9)
        self.T = np.linspace(900, 2000, 12)

    def test_solidus(self):
        tsol = np.array([1130.35520902, 980.6913412, 988.03166758, 1003.60116988,
                         1017.14227837, 1025.892923, 1029.10360679, 1026.73065652,
                         1019.00102304])
        for i in range(len(self.P)):
            self.assertAlmostEqual(tsol[i], self.hlz.TSolidus(self.P[i]))

    def test_liquidus(self):
        tliq = np.array([1885.29107647, 1898.28001541, 1911.26789493, 1924.25471815,
                         1937.24048821, 1950.22520822, 1963.20888128, 1976.19151047,
                         1989.17309886])
        for i in range(len(self.P)):
            self.assertAlmostEqual(tliq[i], self.hlz.TLiquidus(self.P[i]))

    def test_meltFraction(self):
        f = np.array([0., 0.01745679, 0.09695647, 0.17791644, 0.27174548,
                      0.32899751, 0.41441413, 0.52590366, 0.66227606, 0.82198243, 1., 1.])
        for i in range(len(self.T)):
            self.assertAlmostEqual(f[i], self.hlz.F(1.0, self.T[i]))


class test_waterSaturatedSolidus(unittest.TestCase):
    def setUp(self):
        lz = m.lithologies.matthews.klb1()
        self.hlz = m.hydrousLithology(lz, 10.0)
        self.P = np.linspace(0, 8, 9)

    def test_solidus(self):
        tsol = np.array([1130.35520902, 980.6913412, 988.03166758, 1003.60116988,
                         1017.14227837, 1025.892923, 1029.10360679, 1026.73065652,
                         1019.00102304])
        for i in range(len(self.P)):
            self.assertAlmostEqual(tsol[i], self.hlz.TSolidus(self.P[i]))


class test_continuousMelting(unittest.TestCase):
    def setUp(self):
        lz = m.lithologies.matthews.klb1()
        self.hlz = m.hydrousLithology(lz, 1.0, continuous=True)
        self.T = np.linspace(900, 2000, 12)

    def test_meltFraction(self):
        f = np.array([0., 0.01644492, 0.0218464, 0.03096729, 0.06489357, 0.27673482, 0.36147583,
                      0.4837533, 0.63427891, 0.80849053, 1., 1.])
        for i in range(len(self.T)):
            self.assertAlmostEqual(f[i], self.hlz.F(1.0, self.T[i]))


class test_adiabaticMelt(unittest.TestCase):
    def setUp(self):
        lz = m.lithologies.matthews.klb1()
        hlz = m.hydrousLithology(lz, 0.1, continuous=True)
        self.mantle = m.mantle([hlz], [1.0])
        self.col = self.mantle.adiabaticMelt(1350.0)

    def test_adiabaticMelt(self):
        f = np.array([0.00000000e+00, 1.09961923e-05, 2.35590113e-05, 3.66675067e-05,
                      5.01061713e-05, 6.37810468e-05, 7.76402721e-05, 9.16513609e-05,
                      1.05792259e-04, 1.20047132e-04])
        p = np.array([4.33719221, 4.33319221, 4.32919221, 4.32519221, 4.32119221,
                      4.31719221, 4.31319221, 4.30919221, 4.30519221, 4.30119221])
        t = np.array([1437.61516004, 1437.53221586, 1437.4428279, 1437.35316342,
                      1437.26333329, 1437.17338574, 1437.0833474, 1436.99323494,
                      1436.90305967, 1436.81282972])

        for i in range(10):
            self.assertAlmostEqual(f[i], self.col.F.iloc[i])
            self.assertAlmostEqual(p[i], self.col.P.iloc[i])
            self.assertAlmostEqual(t[i], self.col.T.iloc[i])
