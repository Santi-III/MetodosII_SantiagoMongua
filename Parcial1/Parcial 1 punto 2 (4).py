import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
#u funciona como tetha y v como la derivada de tetha respecto al tiempo
def System (r,t,g,L):
    u,w = r
    dudt = w
    dwdt = ((((2*g)/L)-((w*w)*np.cos(u)))*np.sin(u))/((1/3)+(np.sin(u)*np.sin(u)))
    return [dudt,dwdt]

def Integrator2(f,r0,t,g,L):
    h = t[1]-t[0]
    u = np.zeros_like(t)
    w = np.zeros_like(t)

    u[0] = r0[0]
    w[0] = r0[1]

    K1 = np.zeros(2)
    K2 = np.zeros(2)

    for i in range(1,len(t)):
        
        R = np.array([u[i-1],w[i-1]]) #presente
        K1 = f(R,t[i-1],g,L)
        R = np.array([u[i-1]+h*K1[0],w[i-1]+h*K1[1]])
        K2 = f(R,t[i-1]+h,g,L)

        u[i] = u[i-1] + 0.5 * h * (K1[0]+K2[0])
        w[i] = w[i-1] + 0.5 * h * (K1[1]+K2[1])
        
    return u,w

u0 = np.radians(10.0)
r0 = [u0,0.]
g = 9.81
L = 1
N = 100
t = np.linspace(0,0.52,N)

u,w = Integrator2(System,r0, t, g, L)
r = integrate.odeint(System, r0, t, args=(g,L) )
x_runge = (L/2) * np.sin(u)
y_runge = (L/2) * np.cos(u)

x_odeint = (L/2) * np.sin(r[:,0])
y_odeint = (L/2) * np.cos(r[:,0])

plt.figure()
plt.plot(x_runge,y_runge,label = 'Runge 2',color='r')
#plt.plot(x_odeint,y_odeint,label = 'odeint',color='b')
plt.legend()
plt.show()