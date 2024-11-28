import pyMelt as m
from numpy import allclose
import pickle
import os

pyMelt_path = os.path.dirname(os.path.realpath(__file__))[:-6]

f = open(pyMelt_path + '/pyMelt/phaseDiagrams/klb1_holland2018/klb1_holland2018.p', 'rb')
phaseDiagram_object = pickle.load(f)
f.close()
lith = m.phaseDiagramLithology(phaseDiagram_object)

# def test_should_instantiate_phaseDiagramLithology():
#     f = open(pyMelt_path + '/phaseDiagrams/klb1_holland2018/klb1_holland2018.p', 'rb')
#     phaseDiagram_object = pickle.load(f)
#     f.close()
#     lith = m.phaseDiagramLithology(phaseDiagram_object)

def test_should_get_solidusT_from_pdlithology():
    P = 4.0
    Tsolidus = lith.TSolidus(P)
    print(Tsolidus)
    assert allclose(Tsolidus, 1617.9398843848512)

def test_should_get_liquidusT_from_pdlithology():
    P = 4.0
    Tliquidus = lith.TLiquidus(P)
    print(Tliquidus)
    assert allclose(Tliquidus, 1938.6381902319952)

def test_should_get_meltFraction_from_pdlithology():
    P = 4.0
    T = 1650.0
    F = lith.F(P,T)
    print(F)
    assert allclose(F, 0.20553748754164314)

def test_should_return_zero_F_at_subsolidus_T_from_pdlithology():
    P = 4.0
    T = 1100.0
    F = lith.F(P,T)
    print(F)
    assert allclose(F, 0.0)

def test_should_return_unity_F_at_suoperliquidus_T_from_pdlithology():
    P = 4.0
    T = 2100.0
    F = lith.F(P,T)
    print(F)
    assert allclose(F, 1.0)

def test_should_calculate_dTdP_from_pdLithology():
    P = 4.0
    T = 1650.0
    dTdP = lith.dTdP(P, T)
    print(dTdP)
    assert allclose(78.15069929790752, dTdP)

def test_should_calculate_dTdF_from_pdLithology():
    P = 4.0
    T = 1650.0
    dTdF = lith.dTdF(P, T)
    print(dTdF)
    assert allclose(440.1956274636773, dTdF)

def test_should_do_adiabatic_melt_calc_with_pdLithology():
    Tp = 1350.0
    mantle = m.mantle([lith], [1.0], ['lz'])
    column = mantle.adiabaticMelt(Tp)
    print(column.F.iloc[-1])
    # column.plot()
    assert allclose(column.F.iloc[-1], 0.25469077174070326)