"""
==================
Chemistry - Matzen
==================

This module provides the functions implementing the expressions for Mn and Ni partitioning
by Matzen et al. (2017). In most use cases the function `D` is most useful.
"""

import numpy as _np
# Phantom divide by 0 error appears
_np.seterr(divide = 'ignore') 
from pyMelt.core import InputError

def D(state, mineral, element,
      liq_MgO=None, mineral_MgO=None, olv_MgO=None,
      cpx_CaO=None, spn_Cr2O3=None, spn_Al2O3=None, spn_Fe2O3=None,
      **kwargs):
    r"""
    Calculates partition coefficients according to the Matzen (2017) expressions. In default
    operation it will provide partition coefficients for Ni and Mn between liquid and olv, cpx,
    opx, spn, and grt. This function provides an interface to the various KD expressions and
    then calculates D using:

    .. math::
        D_{El}^{min-liq} = \exp \left( \ln \left[ \frac{MgO^{min}}{MgO^{liq}} \right] 
         + \ln K_D^{ol-liq} - \ln K_D^{ol-min} \right)
    
    The constants used in the KD expressions can be provided here and they will propagate
    through the calculation. See the documentation to the KD expressions for the argument
    names.

    Note: At the moment olivine must be in the assemblage. The equations should cancel out
    such that olivine isn't required. Either the equations should be reconfigured so that D
    is calculated directly for mineral-liquid (probably more computationally efficient, but
    at the cost of benchmarking simplicity) or a route through this loop for olivine-absence
    should be created using a fake olivine composition. Simon 19/07/23.

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    mineral : string
        The mineral for which the partition coefficient with melt should be calculated. One
        of 'olv', 'cpx', 'opx', 'spn', or 'grt'.
    element : string
        The element for which the partition coefficient should be calculated. One of 'Ni' 
        or 'Mn'.
    liq_MgO : float, None, or function, default:  None
        The wt% of MgO in the liquid. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'liq_MgO_wtpt'.
    mineral_MgO : float, None, or function, default:  None
        The wt% of MgO in the mineral. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called '<mineral>_MgO_wtpt'. Not required for cpx-Mn.
    olv_MgO : float, None, or function, default:  None
        The wt% of MgO in olivine. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'olv_MgO_wtpt'. Only needed for calculating cpx DMn, it must
        be passed as mineral_MgO when calculating the olivine D.
    cpx_CaO : float, None, or function, default:  None
        The wt% of CaO in cpx. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'cpx_CaO_wtpt'. Only needed for calculating cpx DMn.
    spn_Cr2O3 : float, None, or function, default:  None
        The wt% of Cr2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Cr2O3_wtpt'. Only required for Spn.
    spn_Al2O3 : float, None, or function, default:  None
        The wt% of Al2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Al2O3_wtpt'. Only required for Spn.
    spn_Fe2O3 : float, None, or function, default:  None
        The wt% of Fe2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Fe2O3_wtpt'. Only required for Spn.


    """

    if element == 'Mn':
        kd_olv = kd_olvliq_Mn(state, liq_MgO, **kwargs)
        if mineral == 'olv':
            kd_min = 1.0
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        elif mineral == 'opx':
            kd_min = kd_olvopx_Mn(state, **kwargs)
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        elif mineral == 'cpx':
            Dolv = D_from_kd(1.0, kd_olv, liq_MgO, olv_MgO, state, mineral, **kwargs)
            Dolvcpx = D_olvcpx_Mn(state, cpx_CaO, **kwargs)
            D = Dolv / Dolvcpx
        elif mineral == 'spn':
            kd_min = kd_olvspn_Mn(state, spn_Cr2O3, spn_Al2O3, spn_Fe2O3, **kwargs)
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        elif mineral == 'grt':
            kd_min = kd_olvgrt_Mn(state, **kwargs)
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        else:
            raise InputError("Mineral " + mineral + " not recognised.")
        
    elif element == 'Ni':
        kd_olv = kd_olvliq_Ni(state, **kwargs)
        if mineral == 'olv':
            kd_min = 1.0
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        elif mineral == 'opx':
            kd_min = kd_olvopx_Ni(state, **kwargs)
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        elif mineral == 'cpx':
            kd_min = kd_olvcpx_Ni(state, **kwargs)
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        elif mineral == 'spn':
            kd_min = kd_olvspn_Ni(state, spn_Cr2O3, spn_Al2O3, spn_Fe2O3, **kwargs)
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        elif mineral == 'grt':
            kd_min = kd_olvgrt_Ni(state, **kwargs)
            D = D_from_kd(kd_min, kd_olv, liq_MgO, mineral_MgO, state, mineral, **kwargs)
        else:
            raise InputError("Mineral " + mineral + " not recognised.")
    else:
        raise InputError("Element " + element + " not recognised. It must be one of Ni or Mn "
                         "to be compatible with the Matzen expressions.")

    return D

def D_from_kd(mineral_kd, olv_kd, mineral_MgO=None, liq_MgO=None, state=None, 
              mineral=None, **kwargs):
    r"""
    Converts the KD values and MgO contents of the minerals and melts into the partition
    coefficient.

    .. math::
        D_{El}^{min-liq} = \exp \left( \ln \left[ \frac{MgO^{min}}{MgO^{liq}} \right] 
         + \ln K_D^{ol-liq} - \ln K_D^{ol-min} \right)

    Parameters
    ----------
    mineral_kd : float
        The Kd calculated for element/MgO partitioning between olivine/mineral.
    olv_kd : float
        The Kd calculated for element/MgO partitioning between olivine/liquid
    mineral_MgO : float, None, or function, default:  None
        The wt% of MgO in the mineral. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called '<mineral>_MgO_wtpt'.
    liq_MgO : float, None, or function, default:  None
        The wt% of MgO in the liquid. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'liq_MgO_wtpt'.
    """

    # Extract the liquid MgO content if a function or phase diagram is being used:
    if callable(liq_MgO):
        liq_MgO = liq_MgO(state, **kwargs)
    elif liq_MgO is None and 'phaseDiagram' in kwargs:
        liq_MgO = kwargs['phaseDiagram']('liq_MgO_wtpt',  state)

    # Extract the mineral MgO content if a function or phase diagram is being used:
    if callable(mineral_MgO):
        mineral_MgO = mineral_MgO(state, **kwargs)
    elif mineral_MgO is None and 'phaseDiagram' in kwargs:
        mineral_MgO = kwargs['phaseDiagram'](mineral + '_MgO_wtpt',  state)

    DMg = mineral_MgO / liq_MgO

    return _np.exp(_np.log(DMg) + _np.log(olv_kd) - _np.log(mineral_kd))


def kd_olvliq_Mn(state, liq_MgO=None, matzen_olvliq_Mn_A=0.0088, 
                 matzen_olvliq_Mn_B=-1.503, **kwargs):
    r"""
    Expression for calculating the KD for ol-liq Mn-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Mn-Mg}^{ol-liq} = \exp \left( A MgO^{liq} + B \right)

    where the :math:`MgO^{liq}` is specified in wt%.

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    matzen_olvliq_Mn_A : float, default: 0.0088
        Constant in the KD expression
    matzen_olvliq_Mn_B : float, default: -1.503
        Constant in the KD expression
    liq_MgO : float, None, or function, default:  None
        The wt% of MgO in the liquid. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'liq_MgO_wtpt'.
    """

    if callable(liq_MgO):
        liq_MgO = liq_MgO(state, **kwargs)
    elif liq_MgO is None and 'phaseDiagram' in kwargs:
        liq_MgO = kwargs['phaseDiagram']('liq_MgO_wtpt',  state)
    
    return _np.exp(matzen_olvliq_Mn_A * liq_MgO + matzen_olvliq_Mn_B)



def kd_olvliq_Ni(state, matzen_olvliq_Ni_A=4505.0, 
                 matzen_olvliq_Ni_B=-2.075, **kwargs):
    r"""
    Expression for calculating the KD for ol-liq Ni-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Ni-Mg}^{ol-liq} = \exp \left( \frac{A}{T(K)} + B \right)


    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    matzen_olvliq_Ni_A : float, default: 4505.0
        Constant in the KD expression
    matzen_olvliq_Ni_B : float, default: -2.075
        Constant in the KD expression
    """
    
    return _np.exp(matzen_olvliq_Ni_A / (state['T'] + 273.15) + matzen_olvliq_Ni_B)

def kd_olvopx_Mn(state, matzen_olvopx_Mn_A=-328.0, matzen_olvopx_Mn_B=-0.099,
                 **kwargs):
    r"""
    Expression for calculating the KD for olv-opx Mn-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Mn-Mg}^{olv-opx} = \exp \left( \frac{A}{T(K)} + B \right)


    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    matzen_olvopx_Mn_A : float, default: 328.0
        Constant in the KD expression
    matzen_olvopx_Mn_B : float, default: -0.099
        Constant in the KD expression

    Returns
    -------
    float
        The calculated KD value.
    """
    
    return _np.exp(matzen_olvopx_Mn_A / (state['T'] + 273.15) + matzen_olvopx_Mn_B)

def kd_olvopx_Ni(state, matzen_olvopx_Ni_A=1419, matzen_olvopx_Ni_B=-0.241,
                 **kwargs):
    r"""
    Expression for calculating the KD for olv-opx Ni-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Mn-Mg}^{olv-opx} = \exp \left( \frac{A}{T(K)} + B \right)


    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    matzen_olvopx_Ni_A : float, default: -4171.0
        Constant in the KD expression
    matzen_olvopx_Ni_B : float, default: 1.11
        Constant in the KD expression

    Returns
    -------
    float
        The calculated KD value.
    """
    
    return _np.exp(matzen_olvopx_Ni_A / (state['T'] + 273.15) + matzen_olvopx_Ni_B)

def D_olvcpx_Mn(state, cpx_CaO=None,
                matzen_olvcpx_Mn_A=2.1101e-5, matzen_olvcpx_Mn_B=-7.5872e-4,
                matzen_olvcpx_Mn_C=1.3204e-2, matzen_olvcpx_Mn_D=-0.10348, 
                matzen_olvcpx_Mn_E=1.2147, **kwargs):
    r"""
    Expression for calculating the D for olv-cpx Mn by Matzen et al. (2017) of the
    form:

    .. math::
       D^{olv/cpx}_{Mn} = A (CaO^{cpx})^4 + B (CaO^{cpx})^3 + C (CaO^{cpx})^2 
       + D (CaO^{cpx}) + E


    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    cpx_CaO : float, None, or function, default:  None
        The wt% of CaO in the cpx. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'cpx_CaO_wtpt'.
    matzen_olvcpx_Mn_A : float, default: 2.1101e-5
        Constant in the D expression
    matzen_olvcpx_Mn_B : float, default: -7.5872e-4
        Constant in the D expression
    matzen_olvcpx_Mn_C : float, default: 1.3204e-2
        Constant in the D expression
    matzen_olvcpx_Mn_D : float, default: -0.10348
        Constant in the D expression
    matzen_olvcpx_Mn_E : float, default: 1.2147
        Constant in the D expression

    Returns
    -------
    float
        The calculated D value.
    """

    if callable(cpx_CaO):
        cpx_CaO = cpx_CaO(state, **kwargs)
    elif cpx_CaO is None and 'phaseDiagram' in kwargs:
        cpx_CaO = kwargs['phaseDiagram']('cpx_CaO_wtpt',  state)
    
    D = (matzen_olvcpx_Mn_A * cpx_CaO ** 4 + matzen_olvcpx_Mn_B * cpx_CaO ** 3
         + matzen_olvcpx_Mn_C * cpx_CaO ** 2 + matzen_olvcpx_Mn_D * cpx_CaO
         + matzen_olvcpx_Mn_E) 
    
    return D

def kd_olvcpx_Ni(state, matzen_olvcpx_Ni_A=1773.0, matzen_olvcpx_Ni_B=-0.422,
                 **kwargs):
    r"""
    Expression for calculating the KD for olv-cpx Ni-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Mn-Mg}^{olv-cpx} = \exp \left( \frac{A}{T(K)} + B \right)


    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    matzen_olvcpx_Ni_A : float, default: 1773.0
        Constant in the KD expression
    matzen_olvcpx_Ni_B : float, default: -0.422
        Constant in the KD expression

    Returns
    -------
    float
        The calculated KD value.
    """
    
    return _np.exp(matzen_olvcpx_Ni_A / (state['T'] + 273.15) + matzen_olvcpx_Ni_B)

def kd_olvspn_Mn(state, spn_Cr2O3, spn_Al2O3, spn_Fe2O3, 
                 matzen_olvspn_Mn_A=-2292.2,
                 matzen_olvspn_Mn_B=-0.161, **kwargs):
    r"""
    Expression for calculating the KD for olv-spn Mn-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Mn-Mg}^{olv-spn} = \exp \left( \frac{A Cr\#^{spn}_{wt}}{T(K)} + B \right)

    The :math:`Cr\#` is calculated with:

    .. math::
        Cr\#^{spn}_{wt} = \frac{Cr_2O_3}{Cr_2O_3 + Al_2O_3 + Fe_2O_3}

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    spn_Cr2O3 : float, None, or function, default:  None
        The wt% of Cr2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Cr2O3_wtpt'.
    spn_Al2O3 : float, None, or function, default:  None
        The wt% of Al2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Al2O3_wtpt'.
    spn_Fe2O3 : float, None, or function, default:  None
        The wt% of Fe2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Fe2O3_wtpt'.
    matzen_olvspn_Mn_A : float, default: -2292.2
        Constant in the KD expression
    matzen_olvspn_Mn_B : float, default: -0.161
        Constant in the KD expression

    Returns
    -------
    float
        The calculated KD value.

    """

    if callable(spn_Cr2O3):
        spn_Cr2O3 = spn_Cr2O3(state, **kwargs)
    elif spn_Cr2O3 is None and 'phaseDiagram' in kwargs:
        spn_Cr2O3 = kwargs['phaseDiagram']('spn_Cr2O3_wtpt',  state)

    if callable(spn_Al2O3):
        spn_Al2O3 = spn_Al2O3(state, **kwargs)
    elif spn_Al2O3 is None and 'phaseDiagram' in kwargs:
        spn_Al2O3 = kwargs['phaseDiagram']('spn_Al2O3_wtpt',  state)

    if callable(spn_Fe2O3):
        spn_Fe2O3 = spn_Fe2O3(state, **kwargs)
    elif spn_Fe2O3 is None and 'phaseDiagram' in kwargs:
        spn_Fe2O3 = kwargs['phaseDiagram']('spn_Fe2O3_wtpt',  state)

    # Avoid NAN returns from edges of spinel field in phase diagrams
    if (spn_Al2O3 + spn_Cr2O3 + spn_Fe2O3) < 1e-15:
        return 1.0

    Crn = spn_Cr2O3 / (spn_Cr2O3 + spn_Al2O3 + spn_Fe2O3)

    return _np.exp(matzen_olvspn_Mn_A * Crn / (state['T'] + 273.15) + matzen_olvspn_Mn_B)



def kd_olvspn_Ni(state, spn_Cr2O3, spn_Al2O3, spn_Fe2O3, 
                 matzen_olvspn_Ni_A=1722.0,
                 matzen_olvspn_Ni_B=-1.118, **kwargs):
    r"""
    Expression for calculating the KD for olv-spn Ni-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Ni-Mg}^{olv-spn} = \exp \left( \frac{A Cr\#^{spn}_{wt}}{T(K)} + B \right)

    The :math:`Cr\#` is calculated with:

    .. math::
        Cr\#^{spn}_{wt} = \frac{Cr_2O_3}{Cr_2O_3 + Al_2O_3 + Fe_2O_3}

    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    spn_Cr2O3 : float, None, or function, default:  None
        The wt% of Cr2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Cr2O3_wtpt'.
    spn_Al2O3 : float, None, or function, default:  None
        The wt% of Al2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Al2O3_wtpt'.
    spn_Fe2O3 : float, None, or function, default:  None
        The wt% of Fe2O3 in the spn. If None, it will look for a phaseDiagram object in kwargs
        with a parameter called 'spn_Fe2O3_wtpt'.
    matzen_olvspn_Ni_A : float, default: -2292.2
        Constant in the KD expression
    matzen_olvspn_Ni_B : float, default: -0.161
        Constant in the KD expression

    Returns
    -------
    float
        The calculated KD value.

    """

    if callable(spn_Cr2O3):
        spn_Cr2O3 = spn_Cr2O3(state, **kwargs)
    elif spn_Cr2O3 is None and 'phaseDiagram' in kwargs:
        spn_Cr2O3 = kwargs['phaseDiagram']('spn_Cr2O3_wtpt',  state)

    if callable(spn_Al2O3):
        spn_Al2O3 = spn_Al2O3(state, **kwargs)
    elif spn_Al2O3 is None and 'phaseDiagram' in kwargs:
        spn_Al2O3 = kwargs['phaseDiagram']('spn_Al2O3_wtpt',  state)

    if callable(spn_Fe2O3):
        spn_Fe2O3 = spn_Fe2O3(state, **kwargs)
    elif spn_Fe2O3 is None and 'phaseDiagram' in kwargs:
        spn_Fe2O3 = kwargs['phaseDiagram']('spn_Fe2O3_wtpt',  state)

    # Avoid NAN returns from edges of spinel field in phase diagrams
    if (spn_Al2O3 + spn_Cr2O3 + spn_Fe2O3) < 1e-15:
        return 1.0

    Crn = spn_Cr2O3 / 151.9904 / (spn_Cr2O3 / 151.9904 + spn_Al2O3 / 101.961276 + spn_Fe2O3 / 159.69)

    return _np.exp(matzen_olvspn_Ni_A * Crn / (state['T'] + 273.15) + matzen_olvspn_Ni_B)

def kd_olvgrt_Mn(state, matzen_olvgrt_Mn_A=-4171.0, matzen_olvgrt_Mn_B=1.11,
                 **kwargs):
    r"""
    Expression for calculating the KD for olv-grt Mn-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Mn-Mg}^{olv-grt} = \exp \left( \frac{A}{T(K)} + B \right)


    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    matzen_olvgrt_Mn_A : float, default: -4171.0
        Constant in the KD expression
    matzen_olvgrt_Mn_B : float, default: 1.11
        Constant in the KD expression

    Returns
    -------
    float
        The calculated KD value.
    """
    
    return _np.exp(matzen_olvgrt_Mn_A / (state['T'] + 273.15) + matzen_olvgrt_Mn_B)


def kd_olvgrt_Ni(state, matzen_olvgrt_Ni_A=5264.0, matzen_olvgrt_Ni_B=-2.065,
                 **kwargs):
    r"""
    Expression for calculating the KD for olv-grt Ni-Mg by Matzen et al. (2017) of the
    form:

    .. math::
        K_{D,Ni-Mg}^{olv-grt} = \exp \left( \frac{A}{T(K)} + B \right)


    Parameters
    ----------
    state : pandas.Series
        The state of the system, e.g. temperature (T), pressure (P), melt fraction (F). This
        will most likely be generated automatically by the `MeltingColumn_1D` class.
    matzen_olvgrt_Ni_A : float, default: 5264.0
        Constant in the KD expression
    matzen_olvgrt_Ni_B : float, default: -2.065
        Constant in the KD expression

    Returns
    -------
    float
        The calculated KD value.
    """
    
    return _np.exp(matzen_olvgrt_Ni_A / (state['T'] + 273.15) + matzen_olvgrt_Ni_B)