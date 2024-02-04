import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def f(u,alpha):
    return alpha*u

def Solution(x,alpha):
    return np.exp(alpha*x)

def System(u,t,alpha):
    dudt = alpha*u
    return dudt

def Get(u,t,dt,alpha,u0):
    for k in range(1,len(t)):
        u[k] = ((1+(alpha * t[k])) ** k) * u0
    return u

#Por cuestiones de escala se usaron tiempos 0.11, 0.15 y 0.19 para las gráficas
#La solución es muy inestable por lo cual se toman valores pequeños    
    
t_11 = np.arange(0,2,.11)
t_15 = np.arange(0,2,.15)
t_19 = np.arange(0,2,.19)
dt_11 = t_11[1] - t_11[0]
dt_15 = t_15[1] - t_15[0]
dt_19 = t_19[1] - t_19[0]
u0 = 1
alpha = -1

y_ex = Solution(t_11,alpha)
y_11 = np.zeros_like(t_11)
y_15 = np.zeros_like(t_15)
y_19 = np.zeros_like(t_19)

y_11 = Get(y_11,t_11,dt_11,alpha,u0)
y_15 = Get(y_15,t_15,dt_15,alpha,u0)
y_19 = Get(y_19,t_19,dt_19,alpha,u0)

plt.figure()
plt.plot(t_11,y_ex,label = 'Exacta')
plt.plot(t_11,y_11,label = 't = 1.1 s')
plt.plot(t_15,y_15,label = 't = 1.5 s')
plt.plot(t_19,y_19,label = 't = 1.9 s')
plt.title("EDO soluciones")
plt.xlabel("Tiempo (s)")
plt.ylabel("U")
plt.legend()
plt.show()

