import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_contour(Phi, filename=None, zlabel=r"potential $\Phi$ (V)",
                 cmap=plt.cm.turbo):
    """Plot Phi as a contour plot.
    
    Arguments
    ---------
    Phi : 2D array
          potential on lattice
    filename : string or None, optional (default: None)
          If `None` then show the figure and return the axes object.
          If a string is given (like "contour.png") it will only plot 
          to the filename and close the figure but return the filename.
    cmap : colormap
          pick one from matplotlib.cm          
    """
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(111)

    x = np.arange(Phi.shape[0])
    y = np.arange(Phi.shape[1])
    X, Y = np.meshgrid(x, y)
    Z = Phi[X, Y]
    cset = ax.contourf(X, Y, Z, 20, cmap=cmap)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect(1)

    cb = fig.colorbar(cset, shrink=0.5, aspect=5)
    cb.set_label(zlabel)
    
    if filename:
        fig.savefig(filename)
        plt.close(fig)
        return filename
    else:
        return ax

def plot_surf(Phi, filename=None, offset=-1, zlabel=r'potential $\Phi$ (V)',
             elevation=40, azimuth=-65, cmap=plt.cm.turbo):
    """Plot Phi as a 3D plot with contour plot underneath.
    
    Arguments
    ---------
    Phi : 2D array
          potential on lattice
    filename : string or None, optional (default: None)
          If `None` then show the figure and return the axes object.
          If a string is given (like "contour.png") it will only plot 
          to the filename and close the figure but return the filename.
    offset : float, optional (default: 20)
          position the 2D contour plot by offset along the Z direction
          under the minimum Z value
    zlabel : string, optional
          label for the Z axis and color scale bar
    elevation : float, optional
          choose elevation for initial viewpoint
    azimuth : float, optional
          chooze azimuth angle for initial viewpoint
    cmap : colormap
          pick one from matplotlib.cm
    """
     
    x = np.arange(Phi.shape[0])
    y = np.arange(Phi.shape[1])
    X, Y = np.meshgrid(x, y)
    Z = Phi[X, Y]
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2, linewidth=0.5, color="gray")
    surf = ax.plot_surface(X, Y, Z, cmap=cmap, alpha=0.6)
    cset = ax.contourf(X, Y, Z, 20, zdir='z', offset=offset+Z.min(), cmap=cmap)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel(zlabel)
    ax.set_zlim(offset + Z.min(), Z.max())
    
    ax.view_init(elev=elevation, azim=azimuth)

    cb = fig.colorbar(surf, shrink=0.5, aspect=5)
    cb.set_label(zlabel)
    
    if filename:
        fig.savefig(filename)
        plt.close(fig)
        return filename
    else:
        return ax

def Jacobi_L(Phi):
    """One update in the Jacobi algorithm"""
    Phi[:, 0] = (Phi[:, 1])
    Phi[1:-1, 1:-1] = 0.25*(Phi[2:, 1:-1] + Phi[0:-2, 1:-1] + Phi[1:-1, 2:] + Phi[1:-1, 0:-2])
    return Phi

def Jacobi_P(Phi, rho, dx, dy):
    """One update"""
    Phi[1:-1, 1:-1] = Jacobi_L(Phi)[1:-1, 1:-1]  - 0.25 * dx * dy * rho[1:-1, 1:-1]
    return Phi
    
Nmax = 50
tol = np.finfo(float).eps
max_iter = 100000
u = np.zeros((Nmax, Nmax))
u_old = np.zeros_like(u)

x = np.linspace(0, np.pi, Nmax)
y = np.linspace(0, 1, Nmax)
X, Y = np.meshgrid(x,y)
rho = np.sin(X)*np.cos(np.pi * Y / 2).transpose()

u[0, :] = 0 # x = 0 
u[-1, :] = np.cos(np.pi *y/2) # x = pi
u[:, 0] = 0 # y = 0
u[:, -1] = 0 # y = 1

for k in range(max_iter):
    u_old[:, :] = u
    u = Jacobi_P(u, rho, 1/Nmax, 10/Nmax)
    
    du = np.linalg.norm(u - u_old)
    if du < tol:
        print(f"Converged after {k} iterations")
        break
else:
    print("DNC")


plot_contour(u)
plot_surf(u)

u_exact = (np.sinh((np.pi/2) * X) * np.cos((np.pi/2) * Y))/np.sinh(np.pi**2 /2) - (np.sin(X) * np.cos((np.pi/2) * Y))/(1+np.pi**2/4)
plot_contour(u_exact.transpose())
plot_surf(u_exact.transpose())
plot_contour((u-u_exact.transpose()))
plot_surf((u-u_exact.transpose()))
print(np.linalg.norm(u - u_exact.transpose())/Nmax)
