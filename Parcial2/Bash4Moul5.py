import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def System(y,t):
    dydt = t + 2*y
    return dydt

def IntegratorAdams(f,y0,t):
    
    h = t[1] - t[0]
    
    y = np.zeros_like(t)
    
    y[0] = y0
    
    K11 = 0
    K12 = 0
    K13 = 0
    K14 = 0
    
    #Runge-Kutta 4 orden para los primeros cuatro pasos
    for i in range(4):
        state1 = y[i]
        K11 = f(state1,t[i])
        state2 = y[i]+0.5*h*K11
        K12 = f(state2,t[i]+0.5*h)
        state3 = y[i]+0.5*h*K12
        K13 =  f(state3,t[i]+0.5*h)
        state4 = y[i]+h*K13
        K14 =  f(state4,t[i]+h)
        
        y[i+1] = y[i] + (h/6)*(K11+2*K12+2*K13+K14)
    
    yc = y.copy()
    
    for i in range(5,len(t)):
        
        present = y[i-1]
        past = y[i-2]
        past2 = y[i-3]
        past3 = y[i-4]
        
        K11 = f(present,t[i-1])
        K12 = f(past,t[i-2])
        K13 = f(past2,t[i-3])
        K14 = f(past3,t[i-4])
        K15 = f(past3,t[i-5])
        
        y[i] = y[i-1] + (h/720)*(1901*K11-2774*K12+2616*K13-1274*K14+251*K15)
        
        yc[i] = y[i]
        
        # Futuro
        
        futuro = y[i]
        K16 = f(futuro,t[i])
        
        yc[i] = yc[i-1] + (h/1440)*(475*K16+1427*K15-798*K14+482*K13-173*K12+27*K11)
    
    return y,yc

t = np.linspace(0.,4.,100)
y0 = 0.
sol = odeint(System, y0, t)
mu=1.3
ybash,ymoulton = IntegratorAdams(System,y0,t)
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122)
differences = np.abs(sol[:,0]-ymoulton)

ax.plot(t,ymoulton,label='Adams Moulton')
ax.plot(t,sol,label='Odeint')
ax.legend()
ax1.plot(t,differences,label='Difference')
ax1.loglog()
plt.show()