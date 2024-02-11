import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def System (x,y):
    dydx = (x*y**2) - (2/x) * y - (1 / x**4)
    return dydx

def Integrator2(f,y0,x):
    h = x[2]-x[1]
    y = np.zeros_like(x)
    y[0] = y0

    K1 = 0
    K2 = 0

    for i in range(1,len(x)):
        
        R = np.array([y[i-1]]) #presente
        K1 = f(x[i-1],R)
        R = np.array([y[i-1]+h*K1])
        K2 = f(x[i-1]+h,R)

        y[i] = y[i-1] + 0.5 * h * (K1+K2)
        
    return y

y0 = 0.
x0 = np.sqrt(2)
N = 100
t = np.linspace(1,10,N)
t[0] = x0

y = Integrator2(System, y0, t)
print(t,y)
y_int = integrate.odeint(System, y0, t )
y_real = (t**(-2))

plt.figure()
plt.scatter(t,y,label = 'Runge 2',color='r')
#plt.scatter(t,y_real,label = 'Real')
#plt.plot(t,y_int,label = 'odeint',color='b')
plt.legend()
plt.show()
