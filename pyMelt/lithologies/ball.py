"""
==================
Ball et al. (2022)
==================

Implentation of the modified versions of the Katz et al. (2003) melting models provided by
Ball et al. (2022).
"""

from pyMelt.lithologies.katz import lherzolite as _lherzolite


class primitive_mantle(_lherzolite):
    """
    Implementation of the Ball et al. (2022) modified version of the Katz et al. (2003) model for
    melting primitive mantle

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class.

    - Mcpx:   Mass fraction of cpx in the source. Controls the transition to low-productivity
      harzburgite-type melting.
    - A1:     Parameter used to define solidus.
    - A2:     Parameter used to define solidus.
    - A3:     Parameter used to define solidus.
    - B1:     Parameter used to define lherzolite-liquidus.
    - B2:     Parameter used to define lherzolite-liquidus.
    - B3:     Parameter used to define lherzolite-liquidus.
    - C1:     Parameter used to define liquidus.
    - C2:     Parameter used to define liquidus.
    - C3:     Parameter used to define liquidus.
    - beta1:  Parameter used to calculate melt fraction during cpx-present melting.
    - beta2:  Parameter used to calculate melt fraction during cpx-absent melting.
    - r1:     Parameter used to define cpx reaction coefficient.
    - r2:     Parameter used to define cpx reaction coefficient.

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    also be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: 1187.0
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: 30.0
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt.lithology_class.default_properties['alphaf']
        The thermal expansivity of the melt (1e-6 K-1)
    rhos :       float, default: pyMelt.lithology_class.default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt.lithology_class.default_properties['rhof']
        The density of the melt (kg m-3)
    DeltaS :     float, default: 407.0
        The entropy of fusion J K-1 kg-1
    parameters : dict, default: parameters from Katz et al. (2003)
        The model parameters described above
    """

    def __init__(self,
                 CP=1187.0,
                 alphas=30.0,
                 alphaf=68.0,
                 rhos=3.3,
                 rhof=2.9,
                 DeltaS=407.0,
                 parameters={'Mcpx': 0.1713,
                             'A1': 1085.70,
                             'A2': 132.9,
                             'A3': - 5.1,
                             'B1': 1520.0,
                             'B2': 80.0,
                             'B3': - 3.2,
                             'C1': 1780.0,
                             'C2': 45.0,
                             'C3': - 2.0,
                             'beta1': 1.5,
                             'beta2': 1.2,
                             'r1': 0.9913,
                             'r2': -0.1236
                             }):
        super().__init__(CP=CP,
                         alphas=alphas,
                         alphaf=alphaf,
                         rhos=rhos,
                         rhof=rhof,
                         DeltaS=DeltaS,
                         parameters=parameters)


class depleted_mantle(_lherzolite):
    """
    Implementation of the Ball et al. (2022) modified version of the Katz et al. (2003) model for
    melting depleted mantle

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class.

    - Mcpx:   Mass fraction of cpx in the source. Controls the transition to low-productivity
      harzburgite-type melting.
    - A1:     Parameter used to define solidus.
    - A2:     Parameter used to define solidus.
    - A3:     Parameter used to define solidus.
    - B1:     Parameter used to define lherzolite-liquidus.
    - B2:     Parameter used to define lherzolite-liquidus.
    - B3:     Parameter used to define lherzolite-liquidus.
    - C1:     Parameter used to define liquidus.
    - C2:     Parameter used to define liquidus.
    - C3:     Parameter used to define liquidus.
    - beta1:  Parameter used to calculate melt fraction during cpx-present melting.
    - beta2:  Parameter used to calculate melt fraction during cpx-absent melting.
    - r1:     Parameter used to define cpx reaction coefficient.
    - r2:     Parameter used to define cpx reaction coefficient.

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    also be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: 1187.0
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: 30.0
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt.lithology_class.default_properties['alphaf']
        The thermal expansivity of the melt (1e-6 K-1)
    rhos :       float, default: pyMelt.lithology_class.default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt.lithology_class.default_properties['rhof']
        The density of the melt (kg m-3)
    DeltaS :     float, default: 407.0
        The entropy of fusion J K-1 kg-1
    parameters : dict, default: parameters from Katz et al. (2003)
        The model parameters described above
    """

    def __init__(self,
                 CP=1187.0,
                 alphas=30.0,
                 alphaf=68.0,
                 rhos=3.3,
                 rhof=2.9,
                 DeltaS=407.0,
                 parameters={'Mcpx': 0.1678,
                             'A1': 1085.70,
                             'A2': 132.9,
                             'A3': - 5.1,
                             'B1': 1520.0,
                             'B2': 80.0,
                             'B3': - 3.2,
                             'C1': 1780.0,
                             'C2': 45.0,
                             'C3': - 2.0,
                             'beta1': 1.5,
                             'beta2': 1.2,
                             'r1': 1.2472,
                             'r2': -0.1727
                             }):
        super().__init__(CP=CP,
                         alphas=alphas,
                         alphaf=alphaf,
                         rhos=rhos,
                         rhof=rhof,
                         DeltaS=DeltaS,
                         parameters=parameters)


class mixed_mantle(_lherzolite):
    """
    Implementation of the Ball et al. (2022) modified version of the Katz et al. (2003) model for
    melting 50% depleted and 50% primitive mantle

    To use the same format of parameterisation for another lithology, the parameter values
    may be changed. They are provided as a dictionary during initialisation of the class.

    - Mcpx:   Mass fraction of cpx in the source. Controls the transition to low-productivity
      harzburgite-type melting.
    - A1:     Parameter used to define solidus.
    - A2:     Parameter used to define solidus.
    - A3:     Parameter used to define solidus.
    - B1:     Parameter used to define lherzolite-liquidus.
    - B2:     Parameter used to define lherzolite-liquidus.
    - B3:     Parameter used to define lherzolite-liquidus.
    - C1:     Parameter used to define liquidus.
    - C2:     Parameter used to define liquidus.
    - C3:     Parameter used to define liquidus.
    - beta1:  Parameter used to calculate melt fraction during cpx-present melting.
    - beta2:  Parameter used to calculate melt fraction during cpx-absent melting.
    - r1:     Parameter used to define cpx reaction coefficient.
    - r2:     Parameter used to define cpx reaction coefficient.

    The thermal expansivities, the heat capacity, the densities, and the entropy of fusion may
    also be changed during class initialisation.

    Parameters
    ----------
    CP :         float, default: 1187.0
        The heat capacity (J K-1 kg-1)
    alphas :     float, default: 30.0
        The thermal expansivity of the solid (1e-6 K-1)
    alphaf :     float, default: pyMelt.lithology_class.default_properties['alphaf']
        The thermal expansivity of the melt (1e-6 K-1)
    rhos :       float, default: pyMelt.lithology_class.default_properties['rhos']
        The density of the solid (kg m-3)
    rhof :       float, default: pyMelt.lithology_class.default_properties['rhof']
        The density of the melt (kg m-3)
    DeltaS :     float, default: 407.0
        The entropy of fusion J K-1 kg-1
    parameters : dict, default: parameters from Katz et al. (2003)
        The model parameters described above
    """

    def __init__(self,
                 CP=1187.0,
                 alphas=30.0,
                 alphaf=68.0,
                 rhos=3.3,
                 rhof=2.9,
                 DeltaS=407.0,
                 parameters={'Mcpx': 0.1702,
                             'A1': 1085.70,
                             'A2': 132.9,
                             'A3': - 5.1,
                             'B1': 1520.0,
                             'B2': 80.0,
                             'B3': - 3.2,
                             'C1': 1780.0,
                             'C2': 45.0,
                             'C3': - 2.0,
                             'beta1': 1.5,
                             'beta2': 1.2,
                             'r1': 1.0977,
                             'r2': -0.1437
                             }):
        super().__init__(CP=CP,
                         alphas=alphas,
                         alphaf=alphaf,
                         rhos=rhos,
                         rhof=rhof,
                         DeltaS=DeltaS,
                         parameters=parameters)
