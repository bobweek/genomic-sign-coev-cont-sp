from sympy import besselk, Symbol, oo
from sympy.abc import z,n
import sympy as sp
from numpy import sqrt

x = Symbol('x')

sp.limit(z*besselk(1,z),z,0)

sp.N(besselk(1,sqrt(2)))

sp.N(sp.pi)/10