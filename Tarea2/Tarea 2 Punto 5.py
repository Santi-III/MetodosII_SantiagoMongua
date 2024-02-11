import numpy as np
import matplotlib.pyplot as plt

def System(u,t,q):
    dudt = u**(q)
    return dudt

def Euler3(f,t,q):
    h = t[1]-t[0]
    y = np.zeros_like(t)
    y[0] = 0.1

    K1 = 0
    K2 = 0
    K3 = 0

    for i in range(1,len(t)):
        
        R = np.array([y[i-1]]) #presente
        K1 = f(R,t[i-1],q)
        R = np.array([y[i-1]+h*0.5*K1])
        K2 = f(R,t[i-1]+0.5*h,q)
        R = np.array([y[i-1]-h*K1+2*h*K2])
        K3 = f(R,t[i-1]+h,q)

        y[i] = y[i-1] + (1/6) * h * (K1+4*K2+K3)
        
    return y

N = 100
t = np.linspace(0,1,N)
q = np.array([0.,0.2,0.4,0.7,0.9,1.])

plt.figure()
for i in range (len(q)):
    ye3 = Euler3(System,t,q[i])
    plt.plot(t,ye3,label = 'q = ' + str(q[i]))

plt.legend()
plt.show()