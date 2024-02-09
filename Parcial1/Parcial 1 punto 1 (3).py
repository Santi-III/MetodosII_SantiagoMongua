import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from tqdm import tqdm

class Particle:
    
    def __init__(self, r0, v0, a0, t, m=1, radius=2., g = 9.81, Id=0):
        self.dt = t[1] - t[0]
        # Atributos instantaneos
        self.r = r0
        self.v = v0
        self.a = a0
        
        self.m = m
        self.radius = radius
        self.Id = Id
        self.g = 9.81
        self.p = self.m*self.v
        
        self.f = np.array([0,-self.m*self.g])
        
        self.R = np.zeros((len(t),len(r0)))
        self.V = np.zeros_like(self.R)
        self.A = np.zeros_like(self.R)
        
        self.F = np.zeros_like(self.R)
        
        self.P = np.zeros_like(self.R)
        self.L = np.zeros_like(self.R)
        # Fisica
        self.VEk = np.zeros(len(t))
        self.Ep = 0.
        self.VEp = np.zeros(len(t))
        
    def Evolution(self,i):
        
        self.SetPosition(i)
        self.SetVelocity(i)
        
        self.a = self.f/self.m
        # Euler-Cromer
        self.v += self.dt*self.a
        self.r += self.dt*self.v
    # Setter
    def SetPosition(self,i):
        self.R[i] = self.r
    
    def SetVelocity(self,i):
        self.V[i] = self.v
        self.P[i] = self.m*self.v
        self.F[i] = self.f
        self.L[i] = np.cross(self.r,self.m*self.v)
        self.VEk[i] = 0.5*self.m*np.dot(self.v,self.v)
        Ep = self.m * self.g * self.r[1]
        self.VEp[i] = Ep
    
    # Getter
    def GetPosition(self,scale=1):
        return self.R[::scale]
    
    def GetVelocity(self,scale=1):
        return self.V[::scale]
    
    def GetForce(self,scale=1):
        return self.F[::scale]
 
    def GetMomentum(self,scale=1):
        return self.P[::scale]
    
    def GetKineticEnergy(self,scale=1):
        return self.VEk[::scale]
    
    def GetPotentialEnergy(self,scale=1):
        return self.VEp[::scale]
    
    def GetAngularMomentum(self,scale=1):
        return self.L[::scale]
    
    def CheckLimits(self):
        for i in range(2):
            
            if self.r[i] + self.radius > 20 and self.v[i] > 0.:
                self.v[i] = -(self.v[i]*0.9)
            elif self.r[i] - self.radius < -20 and self.v[i] < 0.:
                self.v[i] = -(self.v[i]*0.9)
                
def GetParticles(N,t):
    
    r0 = np.array([-15.,-10.])
    v0 = np.array([2.,0.])
    a0 = np.array([0.,0.])
    p0 = Particle(r0,v0,a0,t,m=1,radius=2,Id=0)
    return p0

def RunSimulation(t,Particles):
    for it in tqdm(range(len(t)), desc='Running simulation', unit=' Steps' ):
        Particles.Evolution(it)
        Particles.CheckLimits()

    return Particles


dt = 0.001
tmax = 30
t = np.arange(0,tmax,dt)
scale = 100
t1 = t[::scale]
part = GetParticles(1,t)
part = RunSimulation(t,part)

fig = plt.figure(figsize=(40,40))
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

def init():
    
    ax.clear()
    ax.set_xlim(-20,20)
    ax.set_ylim(-20,20)
    
def Update(i):
    
    init()
    ax.set_title(r't =  %.3f s' %(t1[i]))

    KE = 0. # Kinetic energy
     
    x = part.GetPosition(scale)[i,0]
    y = part.GetPosition(scale)[i,1]
        
    vx = part.GetVelocity(scale)[i,0]
    vy = part.GetVelocity(scale)[i,1]
        
    circle = plt.Circle( (x,y), part.radius, color='r', fill=False )
    ax.add_patch(circle)
        
    ax.arrow(x,y,vx,vy,color='k',head_width=0.5,length_includes_head=True)
        
KE = part.GetKineticEnergy(scale)
print(part.GetForce(scale)[:,0])
print(part.GetForce(scale)[:,1])
MomentumT = part.GetMomentum(scale)
EnergyK = part.GetKineticEnergy(scale)
EnergyP = part.GetPotentialEnergy(scale)
Lmomentum = part.GetAngularMomentum(scale)
ax1.scatter(t1, part.GetForce(scale)[:,0],color='k',marker='.',label='Fx')
ax1.scatter(t1, part.GetForce(scale)[:,1],color='r',marker='.',label='Fy')
fig3 = plt.figure(figsize=(10,5))
ax3 = fig3.add_subplot(221)
ax4 = fig3.add_subplot(222)
ax5 = fig3.add_subplot(223)
ax3.plot(t1,MomentumT[:,0],label='px')
ax3.plot(t1,MomentumT[:,1],label='py')
ax4.plot(t1,EnergyK,label='Kinetic')
ax4.plot(t1,EnergyP,label='Potential')
ax4.plot(t1,EnergyK + EnergyP,label='Total')
ax5.plot(t1,Lmomentum,label='Angular')
ax3.legend()
ax4.legend()
Animation = anim.FuncAnimation(fig,Update,frames=len(t1),init_func=init)
plt.show()