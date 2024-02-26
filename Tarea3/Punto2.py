import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint 
import matplotlib.animation as anim
from tqdm import tqdm
#Constants
g = 9.81 
m = 0.1
r = 0.1
d = 0.3
psi_dot = 400

def System(R,t):
    phi,psi,theta,theta_dot = R
    
    I0 = (1/4)*m* r**2 + m * d**2
    Iz = (1/2)*m* r**2
    
    dphi= (Iz*psi_dot*np.cos(y0[2])-Iz*psi_dot*np.cos(theta))/(I0*np.sin(theta)**2 + Iz*np.cos(theta)**2)
    dpsidot = psi_dot
    dthetadot = theta_dot
    dtheta2 = (1/I0)*(dphi**2 * np.sin(theta)*np.cos(theta)*(I0-Iz)-dphi*psi_dot*Iz*np.sin(theta)+m*g*d*np.sin(theta))
    
    return np.array([dphi,dpsidot,dthetadot,dtheta2])

#Condiciones Iniciales
y0 = np.array([0,0,np.pi*(1/4),0])

t = np.linspace(0,8,3000)
scale = 10
solution = odeint(System,y0,t)
phi = solution[:,0]
psi = solution[:,1]
theta = solution[:,2]
theta2 = solution[:,3]

x = d * np.sin(theta) * np.sin(phi)
y = d * np.sin(theta) * np.cos(phi)
z = d * np.cos(theta)

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(projection='3d')

def init():
    
    ax.clear()
    ax.set_xlim(-0.5,0.5)
    ax.set_ylim(-0.5,0.5)
    ax.set_zlim(0,0.5)

def Update(i):
    init()
    ax.set_title(r't =  %.3f segundos' %(t[::scale][i]))
    ax.quiver(0,0,0,x[i],y[i],z[i],color = 'b')
    ax.plot(x[:i],y[:i],z[:i],linestyle='dashed', color='r')
    ax.scatter(x[::scale][i],y[::scale][i],z[::scale][i], color='green')
    
velocities = 0
for i in range(1, len(t)):
    velocities += (phi[i]-phi[i-1]) / (t[i]-t[i-1])
    
vel_average = velocities/len(t)
print("la velocidad de precesión es de: " + str(vel_average))
Animation1 = anim.FuncAnimation(fig,Update,frames=len(t),init_func=init)
plt.show()