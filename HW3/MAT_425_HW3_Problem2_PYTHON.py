import numpy as np
import matplotlib.pyplot as plt

# conditions
x = np.arange(0,5,1)
y = np.array([0,1,4,1,0])
h = 1
n = len(x)-2

## matrix setup ##
A = np.diag(4*np.ones(n)) + np.diag(np.ones(n-1), -1) + np.diag(np.ones(n-1), 1)
A[0,0] = 5
A[-1,-1] = 5

b = np.zeros(n)

for i in range(n):
    b[i] = (6/(h**2)) * (y[i] - 2*y[i+1] + y[i+2])

fxx = np.linalg.solve(A,b)

f_xx = np.zeros(n+2)
f_xx[0] = fxx[0]
f_xx[-1] = fxx[n-1]
f_xx[1:-1] = fxx

xx = np.arange(np.min(x), np.max(x)+0.01, 0.01)
f = np.zeros((n+1,len(xx)))

for i in range(n+1):
    f[i,:] =  (f_xx[i]/6)*((x[i+1]-xx)**3 - (x[i+1]-xx)) + (f_xx[i+1]/6)*((xx-x[i])**3 - (xx-x[i]))+y[i]*(x[i+1]-xx) + y[i+1]*(xx-x[i])

test1 = f[0,0:101]
test2 = xx[0:101]
plt.figure()
plt.title("Cubic Spline Interpolation with Quadratic Run-Out Conditions")
plt.plot(x ,y, 'ko', label = "Data")
plt.plot(xx[0:101], f[0,0:101], label = r'$f_0$')
plt.plot(xx[101:201], f[1,101:201], label = r'$f_1$')
plt.plot(xx[201:301], f[2,201:301], label = r'$f_2$')
plt.plot(xx[301:401], f[3,301:401], label = r'$f_3$')
plt.legend(loc = 'best')
plt.savefig("CubicInterp.png")
