import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib as mpl
import matplotlib.animation as animation
from tqdm import tqdm

Nt = 500
Nx = 60
Ny = 60

x = np.linspace(-5,5,Nx)
y = np.linspace(-5,5,Nx)
t = np.linspace(0,10,Nt)

deltax = x[1] - x[0]
deltay = y[1] - y[0]
deltat = t[1] - t[0]

v = 0.3

def ui(x,y):
    return 5*np.exp(-(x**2 + y**2))

def Initu():
    
    u = np.zeros((Nt,Nx,Ny))
    u[0,:,:] = ui(x,y)
    u[:,:,0] = 0.
    u[:,0,:] = 0.
    u[:,:,-1] = 0.
    u[:,-1,:] = 0.
    return u

u = Initu()

def GetSolution():
    
    for l in tqdm(range(1,len(t))):
        
        #T[l,0,:] = np.sin(20*t[l])
        
        for i in range(1,len(x)-1):
            for j in range(1,len(y)-1):
                u[l,i,j] = (((v*deltat)/(deltax**2)) * (u[l-1,i+1,j] - 2*u[l-1,i,j] + u[l-1,i-1,j])) +\
                ((v*deltat/(deltay**2)) * (u[l-1,i,j+1] - 2*u[l-1,i,j] + u[l-1,i,j-1])) - \
                (((deltat* u[l-1,i,j]) / (2*deltax)) * (u[l-1,i+1,j] - u[l-1,i-1,j])) - \
                (((deltat* u[l-1,i,j]) / (2*deltay)) * (u[l-1,i,j+1] - u[l-1,i,j-1])) + u[l-1,i,j]

GetSolution()

fig = plt.figure(figsize=(15,6))
ax = fig.add_subplot(122, projection='3d')
ax1 = fig.add_subplot(121)
X,Y = np.meshgrid(x,y)

def init():
    ax.set_xlim3d(-5,5)
    ax.set_ylim3d(-5,5)
    ax.set_zlim3d(0,10)

def Update(i):

    ax.clear()
    init()
    
    ax.plot_surface(X,Y,u[i,:,:],cmap='viridis')
    ax1.contourf(X,Y,u[i,:,:].T)
    
Animation = animation.FuncAnimation(fig,Update,frames=len(t),init_func=init)
plt.show()