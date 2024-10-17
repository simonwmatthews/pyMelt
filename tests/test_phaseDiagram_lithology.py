import pyMelt as m
from numpy import allclose
import pickle
import os

pyMelt_path = os.path.dirname(os.path.realpath(__file__))[:-5]

f = open(pyMelt_path + '/phaseDiagrams/klb1_holland2018/klb1_holland2018.p', 'rb')
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
    assert allclose(Tsolidus, 1615.748554732505)

def test_should_get_liquidusT_from_pdlithology():
    P = 4.0
    Tsolidus = lith.TLiquidus(P)
    print(Tsolidus)
    assert allclose(Tsolidus, 1938.6381902319952)