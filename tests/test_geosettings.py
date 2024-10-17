"""
Geosetting Tests

This module provides testing routines for the GeoSettings classes.
"""

import unittest
import pyMelt as m
from numpy import allclose

def test_should_create_an_oib_from_a_supersolidus_start():
    lz = m.lithologies.matthews.klb1()
    mantle = m.mantle([lz], [1.0], ['lz'])
    column = mantle.adiabaticMelt(2000.0, Pstart=8.0)
    column.calculateMajorOxides()
    column.calculateMineralProportions()
    column.calculateTraceElements(c0={'lz':{'La':1.0, 'Yb':1.0}})
    column.calculateStableIsotopes(
        species= 'MgO',
        fractionationFactors= {'olv': 0.1, 
                               'cpx':-0.1, 
                               'opx': 0.05, 
                               'grt':0.0},
        isotopeRatioLabel='d26Mg',
    )

    oib = m.geosettings.intraPlate(column, P_lithosphere=2.0)
    print(oib.chemistry.liq_d26Mg)
    print(oib.chemistry.liq_La)
    assert allclose(oib.chemistry.liq_d26Mg, -0.012978970723878424)
    # The value currently being returned is way too small. There is something
    # going wrong with the homogenization for a supersolidus start
    assert allclose(oib.chemistry.liq_La, 0.0)

def test_should_create_a_mor_from_a_supersolidus_start():
    lz = m.lithologies.matthews.klb1()
    mantle = m.mantle([lz], [1.0], ['lz'])
    column = mantle.adiabaticMelt(2000.0, Pstart=8.0)
    column.calculateMajorOxides()
    column.calculateMineralProportions()
    column.calculateTraceElements(c0={'lz':{'La':1.0, 'Yb':1.0}})
    column.calculateStableIsotopes(
        species= 'MgO',
        fractionationFactors= {'olv': 0.1, 
                               'cpx':-0.1, 
                               'opx': 0.05, 
                               'grt':0.0},
        isotopeRatioLabel='d26Mg',
    )

    mor = m.geosettings.spreadingCentre(column, P_lithosphere=2.0)
    print(mor.chemistry.liq_d26Mg)
    print(mor.chemistry.liq_La)

    assert allclose(mor.chemistry.liq_d26Mg, -0.029178154557631463)
    assert allclose(mor.chemistry.liq_La, 0.0)

def test_should_only_place_liq_comps_in_geosetting():
    lz = m.lithologies.matthews.klb1()
    mantle = m.mantle([lz], [1.0], ['lz'])
    column = mantle.adiabaticMelt(2000.0, Pstart=8.0)
    column.calculateMajorOxides()
    column.calculateMineralProportions()
    column.calculateTraceElements(c0={'lz':{'La':1.0, 'Yb':1.0}})
    column.calculateStableIsotopes(
        species= 'MgO',
        fractionationFactors= {'olv': 0.1, 
                               'cpx':-0.1, 
                               'opx': 0.05, 
                               'grt':0.0},
        isotopeRatioLabel='d26Mg',
    )

    oib = m.geosettings.intraPlate(column, P_lithosphere=2.0)

    for item in oib.chemistry.keys():
        if item[:3] != 'liq':
            print(item)
            assert False
    assert True


# Old tests - but they seem to pass with new pyMelt and pytest

class test_spreadingCentre(unittest.TestCase):
    def setUp(self):
        self.lz = m.lithologies.matthews.klb1()
        self.px = m.lithologies.matthews.kg1()

        self.mantle = m.mantle([self.lz, self.px], [2.0, 1.0], ['lz', 'px'])
        self.column = self.mantle.adiabaticMelt(1350.0)

    def test_default_sc(self):
        setting = m.geosettings.spreadingCentre(self.column)

        self.assertAlmostEqual(setting.tc, 11.161347340791947, places=6,
                               msg="Crustal thickness is incorrect.")

        contributions = {'lz': 0.059825894651427555, 'px': 0.9401741053485725}

        for lith in ['lz', 'px']:
            self.assertAlmostEqual(setting.lithology_contributions[lith], contributions[lith],
                                   places=6, msg="Contribution from " + lith + " is incorrect.")

        self.assertAlmostEqual(setting.P_base_of_crust, 0.3610456279634042, places=6,
                               msg="Pressure at base of crust is incorrect.")

        self.assertAlmostEqual(setting.T.iloc[-1], 1224.1722449731217, places=4,
                               msg="Temperature probably not truncated correctly.")

        self.assertAlmostEqual(setting.F.iloc[-1], 0.273574, places=4,
                               msg="F probably not truncated correctly.")

        self.assertAlmostEqual(setting.lithologies['lz'].F.iloc[-1], 0.072851, places=4,
                               msg="Lz F probably not truncated correctly.")

        self.assertAlmostEqual(setting.lithologies['px'].F.iloc[-1], 0.6749656950797016, places=4,
                               msg="Px F probably not truncated correctly.")

    def test_Tcrys(self):
        setting = m.geosettings.spreadingCentre(self.column)
        Tcryslo, Tcryshi = setting.meltCrystallisationT()

        self.assertAlmostEqual(Tcryslo, 1224.0882859331216, places=6,
                               msg="Tcrys not calculated correctly")

        self.assertAlmostEqual(Tcryshi, 1306.304503376189, places=6,
                               msg="Tcrys not calculated correctly")

    def test_continentalRift(self):
        setting = m.geosettings.spreadingCentre(self.column, P_lithosphere=1.0, extract_melt=True)

        self.assertAlmostEqual(setting.tc, 5.453207375896644, places=6,
                               msg="Crustal thickness is incorrect.")

        self.assertAlmostEqual(setting.P_base_of_crust, 0.9998344279634046, places=6,
                               msg="Pressure at base of crust is incorrect.")

    def test_weighting(self):
        setting = m.geosettings.spreadingCentre(
            self.column, weightingFunction=m.geosettings.weighting_expdecay,
            weighting_wavelength=2.0, weighting_amplitude=0.1)

        self.assertAlmostEqual(setting.tc, 11.744518088013182, places=6,
                               msg="Crustal thickness is incorrect.")

        contributions = {'lz': 0.05714747169325495, 'px': 0.9428525283067452}

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

        contributions = {'lz': 0.028690931456415023, 'px': 0.971309068543585}

        for lith in ['lz', 'px']:
            self.assertAlmostEqual(setting.lithology_contributions[lith], contributions[lith],
                                   places=6, msg="Contribution from " + lith + " is incorrect.")

        self.assertIsNone(setting.melt_flux)

    def test_bouyancy(self):
        setting = m.geosettings.intraPlate(self.column, 1.0, relative_density=0.1)

        self.assertAlmostEqual(setting.melt_flux, 0.6771734157796051, places=6,
                               msg="Melt flux not calculated correctly.")

    def test_weighting(self):
        setting = m.geosettings.intraPlate(self.column, 1.0, relative_density=0.1,
                                           weightingFunction=m.geosettings.weighting_expdecay,
                                           weighting_wavelength=1.0)

        self.assertAlmostEqual(setting.melt_flux, 0.4016489976425273, places=6,
                               msg="Melt flux not calculated correctly.")

        contributions = {'lz': 0.028690931456415023, 'px': 0.971309068543585}

        for lith in ['lz', 'px']:
            self.assertAlmostEqual(setting.lithology_contributions[lith], contributions[lith],
                                   places=6, msg="Contribution from " + lith + " is incorrect.")
