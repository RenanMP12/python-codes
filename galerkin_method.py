# Libraries
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# SYMBOLS
x, y, C1, C2 = symbols('x y C1 C2')

# GALERKIN
y = C1*x*(2 - x) + C2*x*(1 - x)**2 + 1*x

# POLYNOMIAL DERIVATIVES
dy   = y.diff(x)
ddy  = y.diff(x,x)

# RESIDUAL
R = -ddy - y + x**2 
phi1 = y.diff(C1) 
phi2 = y.diff(C2)

# EQUATION SYSTEM
I1 = integrate(R*phi1, (x, 0, 1))
I2 = integrate(R*phi2, (x, 0, 1))
system = [I1, I2]
var = [C1, C2] 
sol = solve(system, var)
C1_ = sol[C1]
C2_ = sol[C2]

# FINAL POLYNOMIAL
y1 = y.subs([(C1,C1_), (C2,C2_)]).simplify()
dy1 = dy.subs([(C1,C1_), (C2,C2_)]).simplify()
ddy1 = ddy.subs([(C1,C1_), (C2,C2_)]).simplify()

# VERIFY BOUNDARY CONDITIONS
print("y(0) = 0: ", y1.subs(x,0).simplify() == 0)
print("dy(1) = 1: ", dy1.subs(x,1).simplify() == 1)

# EXACT SOLUTION
x = Symbol('x')
f = Function('f')(x)
ics = {f.subs(x,0): 0, f.diff(x).subs(x,1): 1}
display(ics)
ode = Eq(-f.diff(x,x) - f + x**2,0)
display(ode)
fsol = dsolve(ode, f, ics = ics)
display(fsol.simplify())

f = x**2 - sin(x)/cos(1) + 2*sin(x)*tan(1) + 2*cos(x) - 2

# PLOT 
fig, ax = plt.subplots(figsize = (12, 7))
x_range = np.arange(0.0, 1.1, 0.05)

x_ = list()
y_galerkin = list()
y_exact = list()

for i in x_range:
    x_.append(i)
    value1 = y1.subs([(x, i)])
    y_galerkin.append(value1)
    
    value2 = f.subs([(x, i)])
    y_exact.append(value2)
    
plt.plot(x_, y_galerkin, '-s')
plt.plot(x_, y_exact, '-*')
plt.xlim([0, 1])
plt.legend(["$Galerkin$", "$Exact \, solution$"])
plt.grid()
plt.show()
fig.savefig('galerkin.png',dpi=300)        # save figure as png (dpi = number of pixels)
