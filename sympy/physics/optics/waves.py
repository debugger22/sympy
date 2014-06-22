"""
This module has all the classes and functions related to waves in optics.

**Contains**

* TWave
* Interference
"""

from __future__ import print_function, division

__all__ = ['TWave', 'Interference']

from sympy import (sympify, pi, sin, cos, sqrt, simplify, Symbol, S, C, I,
    symbols, Derivative)
from sympy.core.expr import Expr
from sympy.physics.units import c
from sympy.geometry import Point, Line


class TWave(Expr):

    r"""
    This is a simple transverse sine wave travelling in a one dimensional space.
    Basic properties are required at the time of creation of the object but
    they can be changed later with respective methods provided.

    It has been represented as :math:`A \times cos(k*x - \omega \times t + \phi )`
    where :math:`A` is amplitude, :math:`\omega` is angular velocity, :math:`k`is
    wavenumber, :math:`x` is a spatial variable to represent the position on the
    dimension on which the wave propagates and :math:`\phi` is phase angle of the wave.


    Arguments
    =========

    amplitude : Sympifyable
        Amplitude of the wave.
    frequency : Sympifyable
        Frequency of the wave.
    phase : Sympifyable
        Phase angle of the wave.
    time_period : Sympifyable
        Time period of the wave.
    n : Sympifyable
        Refractive index of the medium.

    Raises
    =======

    ValueError : When neither frequency nor time period is provided
        or they are not consistent.
    TypeError : When anyting other than TWave objects is added.


    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.physics.optics import TWave
    >>> A1, phi1, A2, phi2, f = symbols('A1, phi1, A2, phi2, f')
    >>> w1 = TWave(A1, f, phi1)
    >>> w2 = TWave(A2, f, phi2)
    >>> w3 = w1 + w2  # Superposition of two waves
    >>> w3
    TWave(sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2), f, phi1 + phi2)
    >>> w3.amplitude
    sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2)
    >>> w3.phase
    phi1 + phi2
    >>> w3.speed
    299792458*m/(n*s)
    >>> w3.angular_velocity
    2*pi*f

    """

    def __init__(
            self,
            amplitude,
            frequency=None,
            phase=S.Zero,
            time_period=None,
            n=Symbol('n')):
        frequency = sympify(frequency)
        amplitude = sympify(amplitude)
        phase = sympify(phase)
        time_period = sympify(time_period)
        n = sympify(n)
        self._frequency = frequency
        self._amplitude = amplitude
        self._phase = phase
        self._time_period = time_period
        self._n = n
        if time_period is not None:
            self._frequency = 1/self._time_period
        if frequency is not None:
            self._time_period = 1/self._frequency
            if time_period is not None:
                if frequency != 1/time_period:
                    raise ValueError("frequency and time_period should be consistent.")
        if frequency is None and time_period is None:
            raise ValueError("Either frequency or time period is needed.")

    @property
    def frequency(self):
        """
        Returns the frequency of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.frequency
        f
        """
        return self._frequency

    @property
    def time_period(self):
        """
        Returns the time period of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.time_period
        1/f
        """
        return self._time_period

    @property
    def wavelength(self):
        """
        Returns wavelength of the wave.
        It depends on the medium of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.wavelength
        299792458*m/(f*n*s)
        """
        return c/(self._frequency*self._n)

    @property
    def amplitude(self):
        """
        Returns the amplitude of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.amplitude
        A
        """
        return self._amplitude

    @property
    def phase(self):
        """
        Returns the phase angle of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.phase
        phi
        """
        return self._phase

    @property
    def speed(self):
        """
        Returns the speed of travelling wave.
        It is medium dependent.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.speed
        299792458*m/(n*s)
        """
        return self.wavelength*self._frequency

    @property
    def angular_velocity(self):
        """
        Returns angular velocity of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.angular_velocity
        2*pi*f
        """
        return 2*pi*self._frequency

    @property
    def wavenumber(self):
        """
        Returns wavenumber of the wave.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave
        >>> A, phi, f = symbols('A, phi, f')
        >>> w = TWave(A, f, phi)
        >>> w.wavenumber
        pi*f*n*s/(149896229*m)
        """
        return 2*pi/self.wavelength

    def __str__(self):
        """String representation of a TWave."""
        from sympy.printing import sstr
        return type(self).__name__ + sstr(self.args)

    __repr__ = __str__

    def __add__(self, other):
        """
        Addition of two waves will result in their superposition.
        The type of interference will depend on their phase angles.
        """
        if isinstance(other, TWave):
            if self._frequency == other._frequency and self.wavelength == other.wavelength:
                return TWave(sqrt(self._amplitude**2 + other._amplitude**2 + 2 *
                                  self.amplitude*other.amplitude*cos(
                                      self._phase - other.phase)),
                             self.frequency,
                             self._phase + other._phase
                             )
        else:
            raise TypeError(type(other).__name__ + " and TWave objects can't be added.")

    def _eval_rewrite_as_sin(self, *args):
        return self._amplitude*sin(self.wavenumber*Symbol('x')
            - self.angular_velocity*Symbol('t') + self._phase + pi/2, evaluate=False)

    def _eval_rewrite_as_cos(self, *args):
        return self._amplitude*cos(self.wavenumber*Symbol('x')
            - self.angular_velocity*Symbol('t') + self._phase)

    def _eval_rewrite_as_pde(self, *args):
        from sympy import Function
        mu, epsilon, x, t = symbols('mu, epsilon, x, t')
        E = Function('E')
        return Derivative(E(x, t), x, 2) + mu*epsilon*Derivative(E(x, t), t, 2)

    def _eval_rewrite_as_exp(self, *args):
        from sympy import C, I
        exp = C.exp
        return self._amplitude*exp(I*(self.wavenumber*Symbol('x')
            - self.angular_velocity*Symbol('t') + self._phase))


class Interference(object):

    """
    Interference of two transverse waves.
    Currently it implements only YDSE(2 slits).

    Assumptions
    ===========

    .. [1] Slit plane is parallel to the screen.

    Parameters
    ==========

    wave : TWave
        A transverse wave.
    S : Point (2D), tuple, list
        Source point
    S1 : Point (2D), tuple, list
        Point corresponding to slit 1
    S2 : Point (2D), tuple, list
        Point corresponding to slit 2
    O : Point, tuple, list (2D)
        Center of the screen
    """

    def __init__(self, wave, S, S1, S2, O):

        if not isinstance(wave, TWave):
            raise TypeError('Inappropriate argument type. wave has to be an instance of TWave.')
        else:
            self.wave = wave

        if not isinstance(S, Point):
            if isinstance(S, tuple) or isinstance(S, list):
                if len(S) != 2:
                    raise ValueError("Only 2D points are supported.")
                self.S = Point(S)
            else:
                raise TypeError("S can be only Point, tuple or list")
        if not isinstance(S1, Point):
            if isinstance(S1, tuple) or isinstance(S1, list):
                if len(S1) != 2:
                    raise ValueError("Only 2D points are supported.")
                self.S1 = Point(S1)
            else:
                raise TypeError("S1 can be only Point, tuple or list")
        if not isinstance(S2, Point):
            if isinstance(S2, tuple) or isinstance(S2, list):
                if len(S2) != 2:
                    raise ValueError("Only 2D points are supported.")
                self.S2 = Point(S2)
            else:
                raise TypeError("S2 can be only Point, tuple or list")
        if not isinstance(O, Point):
            if isinstance(O, tuple) or isinstance(O, list):
                if len(O) != 2:
                    raise ValueError("Only 2D points are supported.")
                self.O = Point(O)
            else:
                raise TypeError("O can be only Point, tuple or list")

            self.slit_line = Line(S1, S2)
            print(self.slit_line)
            print(self.O)
            self.d = self.S1.distance(S2)
            self.D = self.slit_line.distance(self.O)
            self.screen = Line(O, slope=self.slit_line.slope)

    @property
    def fringe_width(self):
        """
        Returns fringe width of the interference pattern.
        This is also the distance between two consecutive
        dark or bright fringes.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave, Interference
        >>> A, f, phi = symbols('A, f, phi')
        >>> w = TWave(A, f, phi)
        >>> Interference(w).fringe_width
        """
        return self.wave.wavelength*self.D/self.d

    def intensity_distribution(self, p):
        """
        Calculates intensity distribution at a point on the screen.
        This assumes that the distances S1P and S2P are very large
        in comparison to the distance S1S2.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.optics import TWave, Interference
        >>> A, f, phi = symbols('A, f, phi')
        >>> w = TWave(A, f, phi)
        >>> model = Interference(w, S=(-4, 0), S1=(-2, 2), S2=(-2, -2), O=(0, 0))
        >>> model.intensity_distribution(Point(0, 4))

        """
        p = Point(point)
        if self.screen.contains(p):
            S_1p = self.S1.distance(p)
            S_2p = self.S2.distance(p)