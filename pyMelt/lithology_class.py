from scipy.misc import derivative

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
