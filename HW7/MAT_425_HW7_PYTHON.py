import scipy.special as spp
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime #Used in runtime calculation

## Functions ##
def G(x, yn, ym1, dt): # G used for Newton-Ralphson Method
    return [x - (4/3)*yn + (1/3)*ym1 - dt*(2/3)*(x**2 - x**3),1 - dt*(2/3)*(2*x - 3*x**2)]

def NR_Solve(y0, yn, ym1, dt, tol, max_iter,j=1): # Newton-Ralphson Method Solver
    y_prev = y0;
    g = G(y_prev, yn, ym1, dt)[0]
    dg = G(y_prev, yn, ym1, dt)[1]
    yp1 = y_prev - g/dg
    iteration = 1

    while np.abs(yp1 - y_prev) > tol:
        y_prev = yp1
        g = G(y_prev, yn, ym1, dt)[0]
        dg = G(y_prev, yn, ym1, dt)[1]
        yp1 = y_prev - g/dg
        iteration += 1
        if iteration > max_iter:
            print(f"Did Not Converge in {max_iter} interations.\n grid value: {j}")
            return None
    return yp1


## Constants ##
for dt in [1,10,100,1000]:
    tic = datetime.now() # starts timer
    d = 0.00001
    t = np.arange(0,(2/d)+dt,dt)
    a = (1/d)-1
    y = np.zeros(len(t))
    y[0] = d
    y[1] = d

    ## BDF2 method to solve y' = y^2 - y^3##
    for i in range(2,len(t)):
        y[i] = NR_Solve(y[(i-1)], y[(i-1)], y[(i-2)], dt, np.finfo(np.float32).eps, 500000,i)
    toc = datetime.now()
    print(f"Took {toc-tic} to solve y' = y^2-y^3, y(0) = {d} using BDF2 with dt = {i}\n")
    y_exact = 1/(spp.lambertw(a*np.exp(a-t)).real+1) # Exact soln used for plotting

    ## Plotting ##
    plt.figure()
    plt.title(r"Solution to $y'=y^2-y^3$ $y(0)=0.00001$")
    plt.plot(t, y, 'r', label = fr'BDF2, $dt$ = {dt}')
    plt.plot(t, y_exact, 'k', label = 'Exact')
    plt.legend(loc = "best")
    plt.savefig(f"MAT_425_BDF2_dt={dt}.png")
