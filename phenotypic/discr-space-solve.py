from sympy.solvers import solve
from sympy import Symbol, oo
import sympy as sp

mH = Symbol('mu_H')
AH = Symbol('A_H')
BH = Symbol('B_H')
wH = Symbol('omega_H')
GH = Symbol('G_H')
NH = Symbol('N_H')
thH = Symbol('theta_H')

m = Symbol('mu')
A = Symbol('A')
B = Symbol('B')
w = Symbol('omega')
G = Symbol('G')
N = Symbol('N')
th = Symbol('theta')

mP = Symbol('mu_P')
AP = Symbol('A_P')
BP = Symbol('B_P')
wP = Symbol('omega_P')
GP = Symbol('G_P')
NP = Symbol('N_P')
thP = Symbol('theta_P')

VH = Symbol('V_H')
VP = Symbol('V_P')
CHP = Symbol('C_HP')

S = sp.MatrixSymbol('Sigma',2,2)
P = sp.MatrixSymbol('Psi',2,2)
L = sp.MatrixSymbol('Lambda',2,2)

r = Symbol('r')
D = Symbol('D')

S = sp.Matrix([[VH,CHP],[CHP,VP]])
P = sp.Matrix([[G*A-wH,G*B],[-G*B,G*A-wP]])
L = sp.Matrix([[G/N,0],[0,G/N]])

Vsol = solve([-2*(GH*AH+wH)*VH-2*GH*BH*CHP+GH/NH,
              -2*(GP*AP+wP)*VP+2*GP*BP*CHP+GP/NP,
              -(GP*AP+wP)*CHP+GP*BP*VH-(GH*AH+wH)*CHP-GH*BH*VP+sp.sqrt(GH*GP/(NH*NP))],[VH,VP,CHP])

geomean = sp.sqrt(sp.simplify(Vsol[VH]*Vsol[VP]))
num = sp.simplify(sp.powdenest(Vsol[CHP], force=True))
rho = num/geomean
rho = sp.powdenest(rho.subs([(GH,G),(GP,G),(AH,A),(AP,A),(NH,N),(NP,N)]),force=True)
n, d = sp.fraction(rho)
sp.simplify(sp.powsimp(d,force=True))
sp.print_latex(sp.collect(n,G**2*N))