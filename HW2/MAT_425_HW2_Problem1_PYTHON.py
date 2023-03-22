import numpy as np
from scipy import optimize

def root_Solve_Bisect(f, a, b, tol, i = 0):
    xm = (a+b)/2
    if np.sign(f(a))==np.sign(f(b)): # No root contained in interval
        print("No root contained, adjusting bounds")
        while f(a)*f(b) > 0:
            a = 2*(a-xm) + xm
            b = 2*(b-xm) + xm
        print(f"New bounds set: [a,b] = [{a},{b}]")

    if np.abs(f(xm)) < tol: # Found root close enough
        print(f"Root found in {i} iterations! x = {xm}")
        return xm
    elif np.sign(f(a)) == np.sign(f(xm)): # adjust left bound
        return root_Solve_Bisect(f, xm, b, tol, i+1)
    elif np.sign(f(b)) == np.sign(f(xm)): # adjust right bound
        return root_Solve_Bisect(f, a, xm, tol, i+1)

f = lambda x: x**3 - 2*x + 2
out = root_Solve_Bisect(f, 3,5, np.finfo(float).eps)
print(f"Error: {out - optimize.fsolve(f, -4)}")
