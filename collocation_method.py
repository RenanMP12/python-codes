# LIBRARIES
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# COLLOCATION METHOD
## SYMBOLS
x, y, C1, C2 = symbols('x y C1 C2')

## THIRD ORDER POLYNOMIAL
y_ = C1*x*(2 - x) + C2*x*(1 - x)**2 + 1*x

## POLYNOMIAL DERIVATIVER
dy_   = y_.diff(x)
ddy_  = y_.diff(x,x)

## RESIDUALS
R = -ddy_ - y_ + x**2 

## EQUATION SYSTEMS
R1 = R.subs(x, 1/5)
R2 = R.subs(x, 4/5)
system = [R1, R2]
var = [C1, C2] 
sol = solve(system, var)
C1_ = sol[C1]
C2_ = sol[C2]

## POLYNOMIAL ACHIEVED BY THE COLLOCATION METHOD
y2 = y_.subs([(C1,C1_), (C2,C2_)]).simplify()
dy2 = dy_.subs([(C1,C1_), (C2,C2_)]).simplify()
ddy2 = ddy_.subs([(C1,C1_), (C2,C2_)]).simplify()
R_c = -ddy2 - y2 + x**2 
display(y2)

## VERIFY BOUNDARY CONDITIONS
print("y(0) = 0: ", y_.subs(x,0).simplify() == 0)
print("dy(1) = 1: ", dy_.subs(x,1).simplify() == 1)

# EXACT SOLUTION
x = Symbol('x')
f = Function('f')(x)
ics = {f.subs(x,0): 0, f.diff(x).subs(x,1): 1}
display(ics)
ode = Eq(-f.diff(x,x) - f + x**2,0)
display(ode)
fsol = dsolve(ode, f, ics = ics)
fsol.simplify()

# POLYNOMIAL - EXACT SOLUTION
f = x**2 - sin(x)/cos(1) + 2*sin(x)*tan(1) + 2*cos(x) - 2

# PLOT 
fig, ax = plt.subplots(figsize = (12, 7))
x_range = np.arange(0.0, 1.1, 0.05)

x_ = list()
y_collocation = list()
y_exact = list()


for i in x_range:
    x_.append(i)
    value_collocation = y2.subs([(x, i)])
    y_collocation.append(value_collocation)
    
    value_exact = f.subs([(x, i)])
    y_exact.append(value_exact)
    
plt.plot(x_, y_exact, '-s')
plt.plot(x_, y_collocation, '-*')
plt.xlim([0, 1])
plt.legend(["$Exact \, Solution$","$Collocation \, method$"])
plt.grid()
plt.show()
fig.savefig('Lista01_Q8a.png',dpi=300)  # save figure as png (dpi = number of pixels)
