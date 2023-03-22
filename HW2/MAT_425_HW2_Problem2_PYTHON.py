import numpy as np
from scipy import optimize

def zero_Solve_Secant(f, x1, x2, tol):
    i = 0
    x3 = None
    while np.abs(x2-x1) > tol:
        i += 1
        x3 = x2 - ((x2-x1)*f(x2))/(f(x2) - f(x1))
        x1 = x2
        x2 = x3
    print(f"Root found in {i} iterations! x = {x2}")
    return x3

f = lambda x: x**3 - 2*x + 2
out = zero_Solve_Secant(f, 3, 5, np.finfo(float).eps)
print(f"Error: {out - optimize.fsolve(f, -4)}")
