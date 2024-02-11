import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return x + 2*y

def Solucion(x):
    return -0.5*x - (1-np.exp(2*x))/4.

def Euler2(f,x,y0):
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

def Euler3(f,x,y0):
    h = x[1]-x[0]
    y = np.zeros_like(x)
    y[0] = y0

    K1 = 0
    K2 = 0
    K3 = 0

    for i in range(1,len(x)):
        
        R = np.array([y[i-1]]) #presente
        K1 = f(x[i-1],R)
        R = np.array([y[i-1]+h*0.5*K1])
        K2 = f(x[i-1]+0.5*h,R)
        R = np.array([y[i-1]-h*K1+2*h*K2])
        K3 = f(x[i-1]+h,R)

        y[i] = y[i-1] + (1/6) * h * (K1+4*K2+K3)
        
    return y

def Euler4(f,x,y0):
    h = x[1]-x[0]
    y = np.zeros_like(x)
    y[0] = y0

    K1 = 0
    K2 = 0
    K3 = 0
    K4 = 0

    for i in range(1,len(x)):
        
        R = np.array([y[i-1]]) #presente
        K1 = f(x[i-1],R)
        R = np.array([y[i-1]+h*0.5*K1])
        K2 = f(x[i-1]+0.5*h,R)
        R = np.array([y[i-1]+0.5*h*K2])
        K3 = f(x[i-1]+0.5*h,R)
        R = np.array([y[i-1]+h*K3])
        K4 = f(x[i-1]+h,R)

        y[i] = y[i-1] + (1/6) * h * (K1+2*K2+2*K3+K4)
        
    return y

N = 10
t = np.linspace(0,1,N)
y0 = 0
yext = Solucion(t)
ye2 = Euler2(f,t,y0)
ye3 = Euler3(f,t,y0)
ye4 = Euler4(f,t,y0)

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(121)
ax.plot(t,yext,label = 'Exacta',color='r')
ax.plot(t,ye2,label = 'Runge 2',color='b')
ax.plot(t,ye3,label = 'Runge 3',color='g')
ax.plot(t,ye4,label = 'Runge 4',color='m')
ax.set_title("Soluciones")
ax.legend()
ax2 = fig.add_subplot(122)
ax2.plot(t,np.abs(ye2-yext),label = 'Runge 2',color='r')
ax2.plot(t,np.abs(ye3-yext),label = 'Runge 3',color='g')
ax2.plot(t,np.abs(ye4-yext),label = 'Runge 4',color='m')
ax2.set_title("Error frente a exacta")
ax2.legend()
plt.show()