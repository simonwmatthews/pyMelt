"""
=================
Lithology classes
=================

The lithology classes module provides the `lithology` class for implementing melting models, and
the `hydrous_lithology` class for converting an anhydrous lithology to a hydrous lithology.
"""
from scipy.misc import derivative
from scipy.optimize import root_scalar
from copy import deepcopy

# Default constant values taken from Katz et al., 2003:
default_properties = {'CP':     1000.0,  # Heat capacity in J Kg-1 K-1
                      'alphas':   40.0,  # Thermal expansivity of the solid (K-1).
                      'alphaf':   68.0,  # Thermal expansivity of the melt (K-1).
                      'rhos':      3.3,  # Density of the solid (g cm-3).
                      'rhof':      2.9,  # Density of the melt (g cm-3).
                      'DeltaS':  300.0,  # Entropy of fusion. (J kg-1 K-1).
                      }


class lithology(object):
    """
    Lithology base class. This class contains all the parameters and methods required to calculate
    the melting behaviour of a single mantle component.

    Attributes
    ----------
    DeltaS:     float
        Entropy of fusion. (J kg-1 K-1). Default is 300.0.
    CP:     float
        Heat capacity (J Kg-1 K-1). Default is 1000.0.
    alphas:     float
        Thermal expansivity of the solid (K-1). Default is 40.0.
    alphaf:     float
        Thermal expansivity of the melt (K-1). Default is 68.0.
    rhos:   float
        Density of the solid (g cm-3). Default is 3.3.
    rhof:   float
        Density of the melt (g cm-3). Default is 2.9.

    """

    def __init__(self,
                 CP=default_properties['CP'],
                 alphas=default_properties['alphas'],
                 alphaf=default_properties['alphaf'],
                 rhos=default_properties['rhos'],
                 rhof=default_properties['rhof'],
                 DeltaS=default_properties['DeltaS'],
                 parameters={}):
        self.CP = CP
        self.alphas = alphas
        self.alphaf = alphaf
        self.rhos = rhos
        self.rhof = rhof
        self.DeltaS = DeltaS
        self.parameters = parameters

    def TSolidus(self, P):
        """
        Returns the temperature of the solidus at any given pressure.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        None
            The default lithology never melts.
        """

        return None

    def TLiquidus(self, P):
        """
        Returns the temperature of the liquidus at any given pressure.

        Parameters
        ----------
        P:  float
            Pressure (GPa).

        Returns
        -------
        None
            The default lithology never melts.
        """
        return None

    def F(self, P, T):
        """
        Returns the melt fraction at any given pressure and temperature.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        float
            Melt fraction.
        """

        return 0.0

    def dTdF(self, P, T, **kwargs):
        """
        Calculates dT/dF(const. P) numerically. For faster and more accurate results redefine
        this method.

        Parameters
        ----------
        P:  float
            Pressure (GPa)
        T:  float
            Temperature (degC)

        Returns
        -------
        float
            dT/dF(const. P) (K).
        """

        def _to_diff(T, P, kwargs={}):
            return self.F(P, T, **kwargs)

        return 1.0/(derivative(_to_diff, T, dx=0.1, args=(P, kwargs)))

    def dTdP(self, P, T, **kwargs):
        """
        Calculates dT/dP(const. F) numerically. For faster and more accurate results redefine
        this method.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        float
            dT/dP(const. F) (K GPa-1).
        """
        dTdF = self.dTdF(P, T, **kwargs)

        def _to_diff(P, T, kwargs={}):
            return self.F(P, T, **kwargs)

        dFdP = derivative(_to_diff, P, args=(T, kwargs))

        return dTdF*dFdP


class hydrous_lithology(object):
    """
    The hydrous lithology class modifies the melting expressions in a given lithology class so
    that water-present melting can be modelled.

    Parameters
    ----------
    lithology : pyMelt.lithology_classes.lithology instance
        The anhydrous lithology for which hydrous melting should be modelled.
    H2O : float
        The water content of the lithology in wt%.
    D : float, default: 0.01
        Partition coefficient for water partitioning between melt and the (bulk) residue.
    K : float, default: 43.0
        Parameter controlling hydrous solidus position, from Katz et al. (2003). See notes.
    gamma : float, default: 0.75
        Parameter controlling the hydrous solidus position, from Katz et al. (2003). Must take a
        value between 0 and 1. See notes.
    chi1 : float, default: 12.0
        Parameter controlling the H2O content of melt at H2O-satuation, in wt% GPa^(-l).
    chi2 : float, default: 1.0
        Parameter controlling the H2O content of melt at H2O-satuation, in wt% GPa^(-1).
    l : float, default: 0.6
        Parameter controlling the H2O content of melt at H2O-satuation.
    continuous : bool, default: False
        Controls whether water extracted to melt is kept or removed from the system. If False,
        melting will be done assuming batch closed-system melting. If True, the water contents
        are modelled according to continuous melting, i.e., near-fractional melting, with the
        porosity controlled by the phi argument.
    phi : float, default: 0.005
        The porosity to use during continuous melting, if applicable. The default value is 0.5%.

    Notes
    -----
    The parameters controlling the position of the hydrous solidus are defined by eqn. 16 from
    Katz et al. (2003):

    .. math::

        \Delta T(X_{H2O}) = K X^\gamma_{H2O}

    Where :math:`X_{H2O}` is the water content in the melt.

    The water content at saturation is controlled by eqn. 17 of Katz et al. (2003):

    .. math::

        X^{sat}_{H2O} = \chi_1 P^\lambda + \chi_2 P


    """

    def __init__(self, lithology, H2O, D=0.01, K=43.0, gamma=0.75, chi1=12.0, chi2=1.0,
                 l=0.6, continuous=False, phi=0.005):
        self.lithology = lithology

        for m in dir(lithology):
            if '__' not in m and m not in ['TSolidus', 'TLiquidus', 'F', 'dTdF', 'dTdP']:
                f = lithology.__getattribute__(m)
                if callable(f):
                    self.__setattr__(m, f.__func__.__get__(self, self.__class__))
                else:
                    self.__setattr__(m, f)

        self.F_function = self.lithology.F.__func__.__get__(self, self.__class__)
        self.H2O = H2O
        self.K = K
        self.gamma = gamma
        self.D = D
        self.chi1 = chi1
        self.chi2 = chi2
        self.l = l
        self.continuous = continuous
        self.phi = phi

    def H2O_saturation(self, P):
        """
        The H2O content of the melt at H2O saturation. This is Eqn 17 of Katz et al. (2003).

        Parameters
        ----------
        P : float
            Pressure in GPa

        Returns
        -------
        float
            The H2O content of the melt at H2O-saturation, in wt%.
        """

        return self.chi1 * P ** self.l + self.chi2 * P

    def melt_H2O(self, F, P=None):
        """
        Returns the H2O content of the melt at a particular melt fraction, assuming batch
        melting, and H2O saturation controlled by the H2O_saturation method. If H2O is
        saturated, it will remain in the bulk system.

        Parameters
        ----------
        F : float
            The melt fraction. Must be between 0 and 1.
        P : float or None, default: None
            The pressure in GPa. Used for checking for H2O saturation. If set to None this
            check will be bypassed.

        Returns
        -------
        float
            The H2O content of the melt, in wt%.
        """

        if self.continuous is False:
            Cl = self.H2O/((1.0 - F)*self.D + F)
        else:
            Cl = (self.H2O / ((1-self.phi)*self.D + self.phi)
                  * (1-F)**((1-self.phi)*(1-self.D)/((1-self.phi)*self.D)))

        # Test for H2O saturation:
        if P is not None:
            H2Osat = self.H2O_saturation(P)
            if H2Osat < Cl:
                Cl = H2Osat

        return Cl

    def TSolidus(self, P, F=0.0):
        """
        Returns the temperature of the solidus at any given pressure. Since the solidus is a
        function of the water content of the bulk, the effective solidus position will shift
        during melting.

        Parameters
        ----------
        P : float
            Pressure (GPa).
        F : float, default: 0.0
            Melt fraction, between 0 and 1.

        Returns
        -------
        Tsol:   float
            Solidus temperature (degC).
        """

        TSolidus = self.lithology.TSolidus(P) - self.K * self.melt_H2O(F, P)**self.gamma

        return TSolidus

    def TLiquidus(self, P, **kwargs):
        """
        Returns the temperature of the liquidus at any given pressure.
        Parameters
        ----------
        P:  float
            Pressure (GPa).
        Returns
        -------
        Tliq:   float
            Liquidus temperature (degC).
        """
        return self.lithology.TLiquidus(P)

    def F(self, P, T, **kwargs):
        """
        Returns the melt fraction of the hydrated lithology, by iteratively finding the value of F
        which provides the amount of water to the melt such that the same value of F is returned
        by the calculation. This is equivalent to Eqn 19 of Katz et al. (2003), but its form will
        depend on the lithology being used.

        The numerical problem is solved using the SciPy method root_scalar from the optimize
        module. The method used is brentq, so that the problem may be constrained between a
        melt fraction of 0 and 1.

        Parameters
        ----------
        P : float
            Pressure in GPa
        T : float
            Temperature in degC.

        Returns
        -------
        float
            The melt fraction.

        """

        F = root_scalar(self._f_to_solve, bracket=[0, 1], args=(P, T))

        if F.flag != 'converged':
            raise pyMelt.core.convergenceError('The melt fraction calculation did not converge.')
            return _np.nan
        else:
            return F.root

    def dTdF(self, P, T, **kwargs):
        """
        Calculates dT/dF(const. P) for the hydrated lithology numerically.

        Parameters
        ----------
        P:  float
            Pressure (GPa)
        T:  float
            Temperature (degC)

        Returns
        -------
        float
            dT/dF(const. P) (K).
        """

        return 1.0/(derivative(self._to_diff_dTdF, T, dx=0.001, args=(P, kwargs)))

    def dTdP(self, P, T, **kwargs):
        """
        Calculates dT/dP(const. F) for the hydrated lithology numerically.

        The method uses the root_scalar method of the scipy.optimize module to find the
        curve of T-P for which F=const. The scipy derivative method is then used to numerically
        differentiate this curve.

        Parameters
        ----------
        P:  float
            Pressure (GPa).
        T:  float
            Temperature (degC).

        Returns
        -------
        float
            dT/dP(const. F) (K GPa-1).
        """
        if self.F(P, T, **kwargs) == 0:
            dTdP = self.alphas/self.rhos/self.CP
        elif self.F(P, T, **kwargs) == 1:
            dTdP = self.alphas/self.rhos/self.CP
        else:
            F = self.F(P, T, **kwargs)

            dTdP = derivative(self._to_diff_dTdP, P, dx=0.001, args=(F, T, kwargs))

        return dTdP

    # The following methods support the numerical differentiation and root-finding methods used
    # above. The reason they are defined as methods, and do not sit within their "parent" methods
    # is to assist with module testing and troubleshooting.

    # This method rearranges the arguments so that F can be numerically differentiated.
    def _to_diff_dTdF(self, T, P, kwargs={}):
        return self.F(P, T, **kwargs)

    # This method is used to find the P-T curve at which F remains constant, for the calculation of
    # dT/dP (at const. F).
    def _hold_constant_F(self, t, P, F, kwargs={}):
        return self.F(P, t, **kwargs) - F

    # This finds the temperature at a given pressure for which the melt fraction is eqal to the
    # value specified. Used for calculating dT/dP (at const. F).
    def _to_diff_dTdP(self, P, F, T, kwargs={}):
        t = root_scalar(self._hold_constant_F, x0=T, x1=T+10, args=(P, F, kwargs)).root
        return t

    # This provides the misfit function when solving for F (Eqn. 19 of Katz et al., 2003)
    def _f_to_solve(self, x, P, T):
        misfit = 0
        fcalc = self.F_function(P, T, F=x)
        misfit = fcalc - x
        return misfit
