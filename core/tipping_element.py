"""tipping_element module

Provides classes for tipping_element objects
"""


class tipping_element:
    """Abstract tipping_element class
    This class provides the interface for tipping_element classes.
    It should not be used but rather inherited from by the concrete 
    tipping_element classes.
    """
    def __init__(self):
        """Constructor"""
        self._type = None
        self._par = {}

        """dx/dt and jacobian diagonal element of tipping element. 
        This method should be overwritten from the 
        concrete tipping_element classes to implement
        the special form of the tipping element.
        """

    def get_type(self):
        return self._type

    def get_par(self):
        return self._par

    def set_par(self, key, val):
        self._par[key] = val

    def dxdt_diag(self):
        """dx/dt diagonal element of tipping element. 
        This method should be overwritten from the 
        concrete tipping_element classes to implement
        the special form of the tipping element.
        """
        return lambda t, x : 0

    def jac_diag(self):
        """jacobian diagonal element of tipping element. 
        This method should be overwritten from the 
        concrete tipping_element classes to implement
        the special form of the tipping element.
        """
        return lambda t, x : 0

    def tip_state(self):
        return lambda x : 0



class cusp(tipping_element):
    """Concrete class for cusp-like tipping element"""
    
    def __init__(self, a, b, c ):
        """Constructor with additional parameters for cusp"""
        super().__init__()
        self._type = 'cusp'
        self._par['a'] = a
        self._par['b'] = b
        self._par['c'] = c

    def dxdt_diag(self):
        """returns callable of dx/dt diagonal element of cusp"""
        return lambda t, x : self._par['a'] * pow(x,3) + self._par['b'] * x \
                                                       + self._par['c']

    def jac_diag(self):
        """returns callable jacobian diagonal element of cusp."""
        return lambda t, x : 3 * self._par['a'] * pow(x,2) + self._par['b']

    def tip_state(self):
        return lambda x : x > 0

class hopf(tipping_element):
    """Concrete class for tipping_elements following the dynamics of a
    Hopfbifurcation.
    Implementation using a representation with polar coordinates:
    dr/dt=(bif_par-r^2)*r*a, dphi/dt=b (time-dependence of angle)"""

    def __init__(self, a, b, c):
        """Constructor with additional parameters for (half a) Hopf
        element"""
        super().__init__()
        self._type = 'hopf'
        self._par['a'] = a
        self._par['b'] = b
        self._par['c'] = c

    def dxdt_diag(self):
        """returns callable of dx/dt diagonal element of Hopf element"""
        return lambda t, r: ( self._par['c'] - pow(r,2) ) * r * self._par['a']

    def jac_diag(self):
        """returns callable jacobian diagonal element of Hopf"""
        return lambda t, r: self._par['a'] * self._par['c'] - \
                            self._par['a'] * 3 * pow(r,2)
    
    def tip_state(self):
        return lambda x: self._par['c'] > 0
