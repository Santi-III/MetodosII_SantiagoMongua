import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
#from tqdm import tqdm 

class Planeta:
    
    def __init__(self, e, a, alpha, t):
        
        self.t = t
        self.dt = t[1] - t[0] # Paso del tiempo
        
        self.e = e # Excentricidad
        self.a_ = a # Semi-eje mayor
        self.alpha = alpha
        self.c = 6.3542E-11
        self.G = 4*np.pi**2 # Unidades gaussianas
        
        self.r = np.zeros(3)
        self.v = np.zeros_like(self.r)
        self.a = np.zeros_like(self.r)
        
        self.r[0] = self.a_*(1-self.e)
        self.v[1] = np.sqrt( self.G*(1+self.e)/(self.a_*(1.-self.e)) )
        
        self.R = np.zeros((len(t),len(self.r)))
        self.V = np.zeros_like(self.R)
        
        # El valor del pasado
        self.rp = self.r
        
    def GetAceleration(self):
        
        d = np.linalg.norm(self.r)
        self.a = -(self.G/d**3)*(1+(self.alpha/d**3))*self.r
        
    def Evolution(self,i):
        
        self.SetPosition(i)
        self.SetVelocity(i)
        self.GetAceleration()
        
        if i==0:
            self.r = self.rp + self.v*self.dt
        else:
            # rp pasado, r presente rf futuro
            self.rf = 2*self.r - self.rp + self.a*self.dt**2
            self.v = (self.rf - self.rp)/(2*self.dt)
            
            self.rp = self.r
            self.r = self.rf
    
    def SetPosition(self,i):
        self.R[i] = self.r
        
    def SetVelocity(self,i):
        self.V[i] = self.v
    
    def GetPosition(self,scale=1):
        return self.R[::scale]
    
    def GetVelocity(self,scale=1):
        return self.V[::scale]
    
    def GetPerihelio(self):
        
        Dist = np.linalg.norm(self.R,axis=1)
        
        timeup = []
        angles = []
        
        for i in range(1,len(Dist)-1):
            if Dist[i] < Dist[i-1] and Dist[i] < Dist[i+1]:
                timeup.append(self.t[i])
                phi = 6 * np.pi * self.G / (self.a_ * (1-(self.e)**2) * (self.c)**2)
                angles.append(phi)
            
        return timeup,angles

def GetPlanetas(t):
    
    Mercurio = Planeta(0.205630,0.387098,1.1E-8,t)
    return Mercurio

dt = 0.001
tmax = 4
t = np.arange(0.,tmax,dt)
Planetas = GetPlanetas(t)


def RunSimulation(t,Planetas):
    
    #for it in tqdm(range(len(t)), desc='Running simulation', unit=' Steps' ):
    for i in range(len(t)):
        Planetas.Evolution(i)
    return Planetas

Planetas = RunSimulation(t,Planetas)
timeup, angles = Planetas.GetPerihelio()
times = []
angle = []
print(angles)
for i in range(len(timeup)):
    if (i%2) != 0:
        times.append(timeup)
        theta = angles[i]*360000
        angle.append(theta)

plt.figure()
plt.plot(times,angle)
plt.show()