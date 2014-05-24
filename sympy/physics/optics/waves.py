from __future__ import print_function, division

__all__ = ['TWave']

from sympy import sympify, pi, cos, Symbol, sqrt, simplify


class TWave(object):

    """
    A transverse wave.


    Examples
    ========

    >>> from sympy import *
    >>> from sympy.physics.optics import *
    >>> A1, f1, phi1, A2, f2, phi2 = symbols('A1, f1, phi1, A2, f2, phi2')
    >>> W1 = TWave('W1', A1, f1, phi1)
    >>> W2 = TWave('W2', A2, f2, phi2)
    """

    def __init__(
            self,
            name,
            amplitude,
            frequency=None,
            phase=None,
            time_period=None):
        if not isinstance(name, str):
            raise TypeError('Invalid name.')
        self._name = name
        self._frequency = sympify(frequency)
        self._amplitude = sympify(amplitude)
        self._phase = sympify(phase)
        self._time_period = sympify(time_period)
        if time_period is not None:
            self._frequency = 1/self._time_period

    @property
    def frequency(self):
        """
        Returns the frequency of the wave.
        """
        return self._frequency

    @property
    def time_period(self):
        """
        Returns the time period of the wave.
        """
        return self._time_period

    @property
    def wavelength(self):
        """
        Returns wavelength of the wave.
        """
        return self._frequency*2*pi

    @property
    def amplitude(self):
        """
        Returns the amplitude of the wave.
        """
        return self._amplitude

    @property
    def phase(self):
        """
        Returns the phase angle of the wave.
        """
        return self._phase

    @property
    def speed(self):
        """
        Returns the speed of travelling wave.
        """
        return self.wavelength*self._frequency

    @property
    def energy(self):
        """
        Returns the energy of the wave.
        """

    @property
    def intensity(self):
        """
        Returns intensity of the wave.
        """

    def __repr__(self):
        return repr(self._amplitude *
                    cos(2 *
                        pi /
                        self._frequency *
                        Symbol('t') +
                        self._phase))

    __str__ = __repr__

    def __add__(self, other):
        return TWave(self._name + other._name,
                     sqrt(self._amplitude**2 + other._amplitude**2 - 2 *
                          self.amplitude*other.amplitude*cos(
                              self._phase - other.phase)),
                     self.frequency,
                     self._phase + other._phase
                     )
