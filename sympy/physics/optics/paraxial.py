"""
This module contains code for paraxial approximation.

**Contains**

* ParaxialRay
"""

from __future__ import division

__all__ = ['ParaxialRay']

from sympy import S, Symbol, sympify, sqrt, Matrix
from sympy.core.basic import Basic
from sympy.geometry.line3d import Ray3D
from .waves import TWave



#################################
#   First order approximation   #
#################################

def _sin(theta):
    return sympify(theta)


def _cos(theta):
    return S.One


def _tan(theta):
    return sympify(theta)

