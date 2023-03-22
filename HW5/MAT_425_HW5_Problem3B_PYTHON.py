import numpy as np
import matplotlib.pyplot as plt

dt = 2.1
a = -1
max_iter = 10
yIE = np.zeros(max_iter+1)
yIE[0] = 1
yCN = np.zeros(max_iter+1)
yCN[0] = 1
t = np.arange(0, 21+dt, dt)


for i in range(max_iter-1):
    yIE[i+1] = yIE[i]/(1-a*dt)
    yCN[i+1] = yCN[i] * ((1 + a*dt/2)/(1 - a*dt/2))

plt.figure()
plt.plot(t, yIE, "r", label = 'I.E.')
plt.plot(t, yCN, "b", label = 'C.N.')
plt.plot(np.linspace(0,20,500), np.exp(-np.linspace(0,21,500)), 'm', label = "Exact")
plt.legend()
plt.savefig(f"IE_CN_dt=2_1")
