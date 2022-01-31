"""
Geosetting Tests

This module provides testing routines for the GeoSettings classes.
"""

import unittest
import pyMelt as m


class test_spreadingCentre(unittest.TestCase):
    def setUp(self):
        self.lz = m.lithologies.matthews.klb1()
        self.px = m.lithologies.matthews.kg1()

        self.mantle = m.mantle([self.lz, self.px], [2.0, 1.0], ['lz', 'px'])
        self.column = self.mantle.adiabaticMelt(1350.0)

    def test_default_sc(self):
        setting = m.geosettings.spreadingCentre(self.column)

        self.assertAlmostEqual(setting.tc, 11.164687478434844, places=6,
                               msg="Crustal thickness is incorrect.")

        contributions = {'lz': 0.05985641132963442, 'px': 0.9401435886703656}

        for lith in ['lz', 'px']:
            self.assertAlmostEqual(setting.lithology_contributions[lith], contributions[lith],
                                   places=6, msg="Contribution from " + lith + " is incorrect.")

        self.assertAlmostEqual(setting.P_base_of_crust, 0.3610456279634042, places=6,
                               msg="Pressure at base of crust is incorrect.")

        self.assertAlmostEqual(setting.T.iloc[-1], 1224.183027, places=4,
                               msg="Temperature probably not truncated correctly.")

        self.assertAlmostEqual(setting.F.iloc[-1], 0.273574, places=4,
                               msg="F probably not truncated correctly.")

        self.assertAlmostEqual(setting.lithologies['lz'].F.iloc[-1], 0.072851, places=4,
                               msg="Lz F probably not truncated correctly.")

        self.assertAlmostEqual(setting.lithologies['px'].F.iloc[-1], 0.675019, places=4,
                               msg="Px F probably not truncated correctly.")

    def test_Tcrys(self):
        setting = m.geosettings.spreadingCentre(self.column)
        Tcryslo, Tcryshi = setting.meltCrystallisationT()

        self.assertAlmostEqual(Tcryslo, 1224.0990681375988, places=6,
                               msg="Tcrys not calculated correctly")

        self.assertAlmostEqual(Tcryshi, 1306.3335694054647, places=6,
                               msg="Tcrys not calculated correctly")

    def test_continentalRift(self):
        setting = m.geosettings.spreadingCentre(self.column, P_lithosphere=1.0, extract_melt=True)

        self.assertAlmostEqual(setting.tc, 5.455312885346313, places=6,
                               msg="Crustal thickness is incorrect.")

        self.assertAlmostEqual(setting.P_base_of_crust, 0.9998344279634046, places=6,
                               msg="Pressure at base of crust is incorrect.")

    def test_weighting(self):
        setting = m.geosettings.spreadingCentre(
            self.column, weightingFunction=m.geosettings.weighting_expdecay,
            weighting_wavelength=1.0, weighting_amplitude=0.001)

        self.assertAlmostEqual(setting.tc, 22.32723842461517, places=6,
                               msg="Crustal thickness is incorrect.")

        contributions = {'lz': 0.017701315986741767, 'px': 0.9822986840132583}

        for lith in ['lz', 'px']:
            self.assertAlmostEqual(setting.lithology_contributions[lith], contributions[lith],
                                   places=6, msg="Contribution from " + lith + " is incorrect.")


class test_intraPlate(unittest.TestCase):
    def setUp(self):
        self.lz = m.lithologies.matthews.klb1()
        self.px = m.lithologies.matthews.kg1()

        self.mantle = m.mantle([self.lz, self.px], [2.0, 1.0], ['lz', 'px'])
        self.column = self.mantle.adiabaticMelt(1350.0)

    def test_default(self):
        setting = m.geosettings.intraPlate(self.column, 1.0)

        self.assertAlmostEqual(setting.F.iloc[-1], 0.175816, places=4,
                               msg="F probably not truncated correctly.")

        contributions = {'lz': 0.028759, 'px': 0.971241}

        for lith in ['lz', 'px']:
            self.assertAlmostEqual(setting.lithology_contributions[lith], contributions[lith],
                                   places=6, msg="Contribution from " + lith + " is incorrect.")

        self.assertIsNone(setting.melt_flux)

    def test_bouyancy(self):
        setting = m.geosettings.intraPlate(self.column, 1.0, relative_density=0.1)

        self.assertAlmostEqual(setting.melt_flux, 0.677310692605866, places=6,
                               msg="Melt flux not calculated correctly.")

    def test_weighting(self):
        setting = m.geosettings.intraPlate(self.column, 1.0, relative_density=0.1,
                                           weightingFunction=m.geosettings.weighting_expdecay,
                                           weighting_wavelength=1.0)

        self.assertAlmostEqual(setting.melt_flux, 0.4017655494238036, places=6,
                               msg="Melt flux not calculated correctly.")

        contributions = {'lz': 0.028759, 'px': 0.971241}

        for lith in ['lz', 'px']:
            self.assertAlmostEqual(setting.lithology_contributions[lith], contributions[lith],
                                   places=6, msg="Contribution from " + lith + " is incorrect.")
