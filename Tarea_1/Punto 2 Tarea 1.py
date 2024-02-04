import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from tqdm import tqdm

class Particle:
    
    def __init__(self, r0, v0, a0, t, m=1, radius=2., Id=0):
        
        self.dt = t[1] - t[0]
        # Atributos instantaneos
        self.r = r0
        self.v = v0
        self.a = a0
        
        self.m = m
        self.radius = radius
        self.Id = Id
        self.p = self.m*self.v
        self.f = self.m*self.a
        # Historial
        
        self.R = np.zeros((len(t),len(r0)))
        self.V = np.zeros_like(self.R)
        self.A = np.zeros_like(self.R)
        
        self.F = np.zeros_like(self.R)
        self.P = np.zeros_like(self.R)
        self.L = np.zeros_like(self.R)
    
        # Fisica
        self.K = 20.
        self.VEk = np.zeros(len(t))
        self.Ep = 0.
        self.VEp = np.zeros(len(t))
        
    def Evolution(self,i):
        
        self.SetPosition(i)
        self.SetVelocity(i)
        
        self.a = self.f/self.m
        
        self.SetPotentialEnergy(i)
        # Euler-Cromer
        self.v += self.dt*self.a
        self.r += self.dt*self.v
        
        
    def CalculateForce(self,p):
        
        d = np.linalg.norm(self.r - p.r)
        
        compresion = self.radius + p.radius - d
        
        if compresion >= 0:
            
            Fn = self.K * compresion**3
            Ep = self.K * (compresion**4) / 4
            self.n = (self.r - p.r)/d     
            self.f = np.add(self.f,Fn*self.n)
            self.Ep += Ep
            # Falta implementar energía potencial 
            
    # Aca debes agregar la energía potencial
    def ResetForce(self):
        self.f[:] = 0.
        self.a[:] = 0.
        self.Ep = 0.
    # Setter
    def SetPosition(self,i):
        self.R[i] = self.r
    
    def SetVelocity(self,i):
        self.V[i] = self.v
        self.P[i] = self.m*self.v
        self.L[i] = np.cross(self.r,self.m*self.v)
        self.VEk[i] = 0.5 * self.m * np.dot(self.v,self.v)
    
    def SetPotentialEnergy(self,i):
        self.VEp[i] = self.Ep
    
    # Getter
    def GetPosition(self,scale=1):
        return self.R[::scale]
    
    def GetVelocity(self,scale=1):
        return self.V[::scale]
 
    def GetMomentum(self,scale=1):
        return self.P[::scale]
    
    def GetKineticEnergy(self,scale=1):
        return self.VEk[::scale]
    
    def GetPotentialEnergy(self,scale=1):
        return self.VEp[::scale]
    
    def GetAngularMomentum(self,scale=1):
        return self.L[::scale]
    # Debes agregar las paredes en este punto

def GetParticles(N,t):
    Particles = []
    for i in range (N):
        r = np.random.uniform(low=-20, high=20, size=2)
        v = np.random.uniform(low=-5, high=5, size=2)
        a = np.array([0.,0.])
        Pi = Particle(r,v,a,t,m=1,radius=2,Id=0)
        Particles.append(Pi)
        #Particles[i] = Pi
    return Particles

def RunSimulation(t,Particles):
    for it in tqdm(range(len(t)), desc='Running simulation', unit=' Steps' ):
        
        for i in range(len(Particles)):
            for j in range(len(Particles)):
                if i!=j:
                    Particles[i].CalculateForce(Particles[j])
        
        for i in range(len(Particles)):
            Particles[i].Evolution(it)
            Particles[i].ResetForce()

    return Particles

dt = 0.001
tmax = 10
t = np.arange(0,tmax,dt)
Particles = GetParticles(10,t)
Particles = RunSimulation(t,Particles)
scale = 200
t1 = t[::scale]
#Particles[0].GetPosition()

fig = plt.figure(figsize=(40,40))
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

def init():
    
    ax.clear()
    ax.set_xlim(-20,20)
    ax.set_ylim(-20,20)
    
def Update(i):
    
    init()
    ax.set_title(r't %.3f' %(t1[i]))
    
    for p in Particles:
        
        x = p.GetPosition()[i,0]
        y = p.GetPosition()[i,1]
        
        vx = p.GetVelocity()[i,0]
        vy = p.GetVelocity()[i,1]
        
        circle = plt.Circle((x,y), p.radius, color='r', fill=False)
        ax.add_patch(circle)
        ax.arrow(x,y,vx,vy,color='k',head_width=0.5)
    #ax.plot(Particles)
       
Animation = anim.FuncAnimation(fig,Update,frames=len(t1),init_func=init)

MomentumT = Particles[0].GetMomentum(scale)
EnergyK = Particles[0].GetKineticEnergy(scale)
EnergyP = Particles[0].GetPotentialEnergy(scale)
Lmomentum = Particles[0].GetAngularMomentum(scale)

EnergyP *= 0.5

for i in range(1,len(Particles)):
    MomentumT = np.add(MomentumT,Particles[i].GetMomentum(scale))
    EnergyK = np.add(EnergyK,Particles[i].GetKineticEnergy(scale))
    EnergyP = np.add(EnergyP,Particles[i].GetPotentialEnergy(scale))
    Lmomentum = np.add(Lmomentum,Particles[i].GetAngularMomentum(scale))

fig3 = plt.figure(figsize=(10,5))
ax3 = fig3.add_subplot(221)
ax4 = fig3.add_subplot(222)
ax5 = fig3.add_subplot(223)
ax3.plot(t1,MomentumT[:,0],label='px')
ax3.plot(t1,MomentumT[:,1],label='py')
ax4.plot(t1,EnergyK,label='Kinetic')
ax4.plot(t1,EnergyP,label='Potential')
ax4.plot(t1,EnergyK + 0.5*EnergyP,label='Total')
ax5.plot(t1,Lmomentum,label='Angular')
ax3.legend()
ax4.legend()

plt.show()
