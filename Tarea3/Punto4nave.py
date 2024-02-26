import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import matplotlib.animation as anim
from tqdm import tqdm

#Constants
G = 6.67e-11
mt = 5.9736e24
rt = 6.3781e6
ml = 0.07349e24
rl = 1.7374e6
d = 3.844e8
omega = 2.6617e-6

#System
def System(R,t):
    r,pr,phi,pphi = R
    rprime = np.sqrt(1+r**2-2*r*np.cos(phi-omega*t))
    mu = ml / mt
    Delta = G * mt / d**3
    drdt = pr
    dprdt = ((pphi**2) / (r**3)) - Delta * ((1/r**2)+((mu/(rprime**3))*(r - np.cos(phi - omega*t))))
    dphidt = (pphi/(r**2))
    dpphidt = - ((Delta * r) / (rprime**3)) * np.sin(phi - omega*t)
    return np.array([drdt,dprdt,dphidt,dpphidt])
#Método Runge-Kutta 4
def RK(f,R,t):
    r0,pr0,phi0,pphi0 = R
    h = t[1] - t[0]
    r = np.zeros_like(t)
    phi = np.zeros_like(t)
    pr = np.zeros_like(t)
    pphi = np.zeros_like(t)
    r[0] = r0
    pr[0] = pr0
    phi[0] = phi0
    pphi[0] = pphi0
    
    K1 = np.zeros(4)
    K2 = np.zeros(4)
    K3 = np.zeros(4)
    K4 = np.zeros(4)
    
    for i in tqdm(range(1, len(t)), desc='Running simulation', unit=' Steps'):
        R = np.array([r[i-1],pr[i-1],phi[i-1],pphi[i-1]])
        K1 = f(R,t[i-1])
        R = np.array([r[i-1]+h*0.5*K1[0],pr[i-1]+h*0.5*K1[1],phi[i-1]+h*0.5*K1[2],pphi[i-1]+h*0.5*K1[3]])
        K2 = f(R, t[i-1]+0.5*h)
        R = np.array([r[i-1]+h*0.5*K2[0],pr[i-1]+h*0.5*K2[1],phi[i-1]+h*0.5*K2[2],pphi[i-1]+h*0.5*K2[3]])
        K3 = f(R, t[i-1]+0.5*h)
        R = np.array([r[i-1]+h*K3[0],pr[i-1]+h*K3[1],phi[i-1]+h*K3[2],pphi[i-1]+h*K3[3]])
        K4 = f(R, t[i-1]+h)

        r[i] = r[i-1] + (1/6) * h * (K1[0]+2*K2[0]+2*K3[0]+K4[0])
        pr[i] = pr[i-1] + (1/6) * h * (K1[1]+2*K2[1]+2*K3[1]+K4[1])
        phi[i] = phi[i-1] + (1/6) * h * (K1[2]+2*K2[2]+2*K3[2]+K4[2])
        pphi[i] = pphi[i-1] + (1/6) * h * (K1[3]+2*K2[3]+2*K3[3]+K4[3])
        
    return r,pr,phi,pphi

#Condiciones iniciales 
v0 = 1.15e4/d
phi0 = 360 * np.pi/180
the = 20 * np.pi/180
r0 = 6.378e6/d

R0 = [r0, v0*np.cos(the-phi0), phi0, v0*r0*np.sin(the-phi0)]

t = np.arange(0,2.8e6)
scale = 1000

#Solución y conversión
sol = RK(System,R0,t)
r_sol = sol[0]
pr_sol = sol[1]
phi_sol = sol[2]
pphi_sol = sol[3]

xl = d*np.cos(omega*t)
yl = d*np.sin(omega*t)
xt = r_sol*d*np.cos(phi_sol)
yt = r_sol*d*np.sin(phi_sol)

#Animation
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot()

time = t[::scale]*(1/(24*3600))

def init():
    ax.clear()
    ax.set_xlim(-6e8,6e8)
    ax.set_ylim(-6e8,6e8)
    
def Update(i):
    
    init()
    ax.set_title(r't =  %.3f días de la tierra' %(time[i]))
    ax.scatter(0,0,color='r')
    ax.scatter(xl[::scale][i],yl[::scale][i],color = 'blue',linestyle = 'dashed')
    ax.scatter(xt[::scale][i],yt[::scale][i],color = 'green',linestyle = 'dashed')
        
Animation = anim.FuncAnimation(fig,Update,frames=len(time),init_func=init)
plt.show()