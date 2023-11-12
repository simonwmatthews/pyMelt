"""
==================
Chemistry - Data
==================

This module provides some useful data for calculating magma compositions.
"""

from dataclasses import dataclass
import pandas as _pd


@dataclass
class oxide_masses:
    """
    The molecular mass of the oxides.
    """
    SiO2: float = 28.085 + 15.999 * 2
    MgO: float = 24.305 + 15.999
    FeO: float = 55.845 + 15.999
    CaO: float = 40.078 + 15.999
    Al2O3: float = 26.982 * 2 + 15.999 * 3
    Na2O: float = 22.99 * 2 + 15.999
    K2O: float = 39.098 * 2 + 15.999
    MnO: float = 54.938 + 15.999
    TiO2: float = 79.867
    P2O5: float = 2 * 30.974 + 5 * 15.999
    Cr2O3: float = 151.992
    NiO: float = 58.693 + 16
    CoO: float = 44.01
    Fe2O3: float = 55.845 * 2 + 15.999 * 3
    H2O: float = 18.02
    CO2: float = 44.01
    O: float = 15.999


@dataclass
class workman05_dmm:
    """
    The trace element concentrations in the depleted MORB mantle from Workman & Hart (2005). All
    concentrations are in ppmw.
    """
    Rb: float = 0.05
    Ba: float = 0.563
    Th: float = 0.0079
    U: float = 0.0032
    Nb: float = 0.1485
    Ta: float = 0.0096
    La: float = 0.192
    Ce: float = 0.550
    Pb: float = 0.018
    Pr: float = 0.107
    Nd: float = 0.581
    Sr: float = 7.664
    Zr: float = 5.082
    Hf: float = 0.157
    Sm: float = 0.239
    Eu: float = 0.096
    Ti: float = 716.3
    Gd: float = 0.358
    Tb: float = 0.070
    Dy: float = 0.505
    Ho: float = 0.115
    Y: float = 3.328
    Er: float = 0.348
    Yb: float = 0.365
    Lu: float = 0.058


@dataclass
class workman05_D:
    """
    The bulk partition coefficients for MORB production from Workman & Hart (2005).
    """
    Rb: float = 1e-5
    Ba: float = 0.00012
    Th: float = 0.001
    U: float = 0.0011
    Nb: float = 0.0034
    Ta: float = 0.0034
    La: float = 0.01
    Ce: float = 0.022
    Pb: float = 0.014
    Pr: float = 0.027
    Nd: float = 0.031
    Sr: float = 0.025
    Zr: float = 0.033
    Hf: float = 0.035
    Sm: float = 0.045
    Eu: float = 0.050
    Ti: float = 0.058
    Gd: float = 0.056
    Tb: float = 0.068
    Dy: float = 0.079
    Ho: float = 0.084
    Y: float = 0.088
    Er: float = 0.097
    Yb: float = 0.115
    Lu: float = 0.120


@dataclass
class stracke03_bsic:
    """
    The trace element concentrations (ppmw) in bulk subducted igneous crust from Stracke et al. (2003).
    """
    Rb: float = 0.57
    Ba: float = 6.59
    Th: float = 0.088
    U: float = 0.027
    Nb: float = 1.95
    Ta: float = 0.124
    La: float = 1.68
    Ce: float = 5.89
    Pb: float = 0.09
    Nd: float = 7.45
    Sr: float = 81.0
    Zr: float = 64.0
    Hf: float = 1.78
    Sm: float = 2.69
    Eu: float = 1.04
    Ti: float = 7735.0
    Gd: float = 4.03
    Dy: float = 5.01
    Y: float = 28.5
    Er: float = 3.13
    Yb: float = 2.99
    Lu: float = 0.45


@dataclass
class palme13_pm:
    """
    The composition of the primitive mantle (ppmw) from Palme & O'Neill (2013).
    """
    Rb: float = 0.605
    Ba: float = 6.85
    Th: float = 0.0849
    U: float = 0.0229
    Nb: float = 0.595
    Ta: float = 0.043
    La: float = 0.6832
    Ce: float = 1.7529
    Pb: float = 0.185
    Pr: float = 0.2657
    Nd: float = 1.341
    Sr: float = 22.0
    Zr: float = 10.3
    Hf: float = 0.3014
    Sm: float = 0.4347
    Eu: float = 0.1665
    Ti: float = 1265.0
    Gd: float = 0.5855
    Tb: float = 0.1075
    Dy: float = 0.7239
    Ho: float = 0.1597
    Y: float = 4.13
    Er: float = 0.4684
    Yb: float = 0.4774
    Lu: float = 0.07083


@dataclass
class palme13_ci:
    """
    Trace element concentrations in a CI chondrite (ppmw) from Palme & O'Neill (2013).
    """
    Rb: float = 2.32
    Ba: float = 2.42
    Th: float = 0.03
    U: float = 0.00810
    Nb: float = 0.283
    Ta: float = 0.015
    La: float = 0.2414
    Ce: float = 0.6194
    Pb: float = 2.62
    Pr: float = 0.09390
    Nd: float = 0.4737
    Sr: float = 7.79
    Zr: float = 3.63
    Hf: float = 0.1065
    Sm: float = 0.1536
    Eu: float = 0.05883
    Ti: float = 447.0
    Gd: float = 0.2069
    Tb: float = 0.03797
    Dy: float = 0.2558
    Ho: float = 0.05644
    Y: float = 1.46
    Er: float = 0.1655
    Yb: float = 0.1687
    Lu: float = 0.02503



@dataclass
class olv_D:
    """
    Trace element partition coefficients between olivine and melt, compiled by Gibson & Geist (2010).
    """
    Rb: float = 0.0003
    Ba: float = 0.000005
    Th: float = 0.00005
    U: float = 0.00038
    Nb: float = 0.0005
    Ta: float = 0.0005
    La: float = 0.0005
    Ce: float = 0.0005
    Pb: float = 0.003
    Pr: float = 0.0008
    Nd: float = 0.00042
    Sr: float = 0.00004
    Zr: float = 0.0033
    Hf: float = 0.0022
    Sm: float = 0.0011
    Eu: float = 0.0016
    Ti: float = 0.015
    Gd: float = 0.0011
    Tb: float = 0.0015
    Dy: float = 0.0027
    Ho: float = 0.0016
    Y: float = 0.0099
    Er: float = 0.013
    Yb: float = 0.020
    Lu: float = 0.020



@dataclass
class opx_D:
    """
    Trace element partition coefficients between orthopyroxene and melt, compiled by Gibson & Geist
    (2010).
    """
    Rb: float = 0.0002
    Ba: float = 0.000006
    Th: float = 0.002
    U: float = 0.002
    Nb: float = 0.004
    Ta: float = 0.004
    La: float = 0.0031
    Ce: float = 0.0040
    Pb: float = 0.009
    Pr: float = 0.0048
    Nd: float = 0.01200
    Sr: float = 0.0007
    Zr: float = 0.013
    Hf: float = 0.03
    Sm: float = 0.0200
    Eu: float = 0.0130
    Ti: float = 0.086
    Gd: float = 0.0130
    Tb: float = 0.0190
    Dy: float = 0.0110
    Ho: float = 0.0065
    Y: float = 0.052
    Er: float = 0.045
    Yb: float = 0.080
    Lu: float = 0.120



@dataclass
class cpx_D:
    """
    Trace element partition coefficients between clinopyroxene and melt, compiled by Gibson & Geist
    (2010).
    """
    Rb: float =  0.0004
    Ba: float =  0.0004
    Th: float =  0.0059
    U: float =  0.0094
    Nb: float =  0.015
    Ta: float =  0.015
    La: float =  0.0490
    Ce: float =  0.0800
    Pb: float =  0.012
    Pr: float =  0.126
    Nd: float =  0.17800
    Sr: float =  0.091
    Zr: float =  0.119
    Hf: float =  0.284
    Sm: float =  0.2930
    Eu: float =  0.3350
    Ti: float =  0.350
    Gd: float =  0.3500
    Tb: float =  0.4030
    Dy: float =  0.4000
    Ho: float =  0.4270
    Y: float =  0.426
    Er: float =  0.420
    Yb: float =  0.400
    Lu: float =  0.376



@dataclass
class grt_D:
    """
    Trace element partition coefficients between garnet and melt, compiled by Gibson & Geist (2010).
    """
    Rb: float = 0.0002
    Ba: float = 0.00007
    Th: float = 0.009
    U: float = 0.028
    Nb: float = 0.015
    Ta: float = 0.015
    La: float = 0.0010
    Ce: float = 0.0050
    Pb: float = 0.005
    Pr: float = 0.014
    Nd: float = 0.05200
    Sr: float = 0.0007
    Zr: float = 0.270
    Hf: float = 0.400
    Sm: float = 0.2500
    Eu: float = 0.4960
    Ti: float = 0.600
    Gd: float = 0.84800
    Tb: float = 1.4770
    Dy: float = 2.2000
    Ho: float = 3.3150
    Y: float = 3.100
    Er: float = 4.400
    Yb: float = 6.600
    Lu: float = 7.100



@dataclass
class spn_D:
    """
    Trace element partition coefficients between spinel and melt, taken from alphaMELTS.
    """
    Rb: float = 0.0001
    Ba: float = 0.0001
    Th: float = 0.0
    U: float = 0.0
    Nb: float = 0.0
    Ta: float = 0.0
    La: float = 0.0100
    Ce: float = 0.0100
    Pb: float = 0.0
    Pr: float = 0.01
    Nd: float = 0.0100
    Sr: float = 0.0
    Zr: float = 0.0
    Hf: float = 0.0
    Sm: float = 0.0100
    Eu: float = 0.0100
    Ti: float = 0.15
    Gd: float = 0.0100
    Tb: float = 0.0100
    Dy: float = 0.0100
    Ho: float = 0.0100
    Y: float = 0.01
    Er: float = 0.0100
    Yb: float = 0.0100
    Lu: float = 0.0100


@dataclass
class plg_D:
    """
    Trace element partition coefficients between plagioclase and melt, compiled by Gibson & Geist
    (2010).
    """
    Rb: float = 0.03
    Ba: float = 0.33
    Th: float = 0.05
    U: float = 0.11
    Nb: float = 0.01
    Ta: float = 0.0
    La: float = 0.2700
    Ce: float = 0.200
    Pb: float = 0.36
    Pr: float = 0.17
    Nd: float = 0.1400
    Sr: float = 2.0
    Zr: float = 0.01
    Hf: float = 0.01
    Sm: float = 0.1100
    Eu: float = 0.7300
    Ti: float = 0.04
    Gd: float = 0.0660
    Tb: float = 0.0600
    Dy: float = 0.0550
    Ho: float = 0.0480
    Y: float = 0.03
    Er: float = 0.0100
    Yb: float = 0.031
    Lu: float = 0.0250


@dataclass
class klb1_MineralProportions:
    """
    Mineral proportions (Wt%) for KLB1 in the garnet-, spinel-, and plagioclase-field (estimated
    from Jennings and Holland, 2015).
    """
    df: _pd.DataFrame = _pd.DataFrame(
            [
                [0.609, 0.125, 0.119, 0.147, 0.000, 0.000],
                [0.597, 0.233, 0.158, 0.000, 0.012, 0.000],
                [0.646, 0.208, 0.076, 0.000, 0.000, 0.070],
            ],
            columns=["olv", "opx", "cpx", "grt", "spn", "plg"],
            index=["grt_field", "spn_field", "plg_field"],
        )


@dataclass
class kg1_MineralProportions:
    """
    Mineral proportions (Wt%) for KG1 in the garnet-, spinel-, and plagioclase-field (estimated
    from Matthews et al., 2021).
    """
    df: _pd.DataFrame = _pd.DataFrame(
            [
                [0.181, 0.012, 0.422, 0.385, 0.000, 0.000],
                [0.110, 0.178, 0.641, 0.000, 0.071, 0.000],
                [0.118, 0.150, 0.655, 0.000, 0.000, 0.067],
            ],
            columns=["olv", "opx", "cpx", "grt", "spn", "plg"],
            index=["grt_field", "spn_field", "plg_field"],
        )


@dataclass
class mo91_MineralProportions:
    """
    Mineral proportions (Wt%) for lherzolite in the garnet-, spinel-, and plagioclase-field, from
    McKenzie & O'Nions (1991).
    """
    df : _pd.DataFrame = _pd.DataFrame(
            [
                [0.598, 0.221, 0.076, 0.115, 0.000, 0.000],
                [0.578, 0.270, 0.119, 0.000, 0.033, 0.000],
                [0.636, 0.263, 0.012, 0.000, 0.000, 0.089],
            ],
            columns=["olv", "opx", "cpx", "grt", "spn", "plg"],
            index=["grt_field", "spn_field", "plg_field"],
        )


@dataclass
class ionic_radii:
    """
    The ionic radii (in m) of the trace elements, as compiled in the alphaMELTS input files.
    """
    Na: float =  1.18e-10
    P: float =  0.38e-10
    K: float =  1.51e-10
    Sc: float =  0.87e-10
    Ti: float =  0.74e-10
    V: float =  0.54e-10
    Cr: float =  0.615e-10
    Mn: float =  0.96e-10
    Co: float =  0.9e-10
    Ni: float =  0.69e-10
    Ga: float =  0.62e-10
    Rb: float =  1.61e-10
    Sr: float =  1.26e-10
    Y: float =  1.019e-10
    Nb: float =  0.64e-10
    Cd: float =  1.10e-10
    In: float =  0.92e-10
    Cs: float =  1.74e-10
    Ba: float =  1.42e-10
    La: float =  1.16e-10
    Ce: float =  1.143e-10
    Pr: float =  1.126e-10
    Nd: float =  1.109e-10
    Sm: float =  1.079e-10
    Eu: float =  1.066e-10
    Gd: float =  1.053e-10
    Tb: float =  1.040e-10
    Dy: float =  1.027e-10
    Ho: float =  1.015e-10
    Er: float =  1.004e-10
    Tm: float =  0.994e-10
    Yb: float =  0.985e-10
    Lu: float =  0.977e-10
    Ta: float =  0.64e-10
    Ra: float =  1.48e-10
    Ac: float =  1.25e-10
    Th: float =  1.035e-10
    Pa: float =  0.78e-10
    U: float =  0.975e-10


@dataclass
class ionic_charge:
    """
    The charges on the trace element ions, as compiled in the alphaMELTS input files.
    """
    Na: int = 1
    P: int = 5
    K: int = 1
    Sc: int = 3
    Ti: int = 4
    V: int = 5
    Cr: int = 3
    Mn: int = 2
    Co: int = 2
    Ni: int = 2
    Ga: int = 3
    Rb: int = 1
    Sr: int = 2
    Y: int = 3
    Nb: int = 5
    Cd: int = 2
    In: int = 3
    Cs: int = 1
    Ba: int = 2
    La: int = 3
    Ce: int = 3
    Pr: int = 3
    Nd: int = 3
    Sm: int = 3
    Eu: int = 3
    Gd: int = 3
    Tb: int = 3
    Dy: int = 3
    Ho: int = 3
    Er: int = 3
    Tm: int = 3
    Yb: int = 3
    Lu: int = 3
    Ta: int = 5
    Ra: int = 2
    Ac: int = 3
    Th: int = 4
    Pa: int = 5
    U: int = 4