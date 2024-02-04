import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from tqdm import tqdm
from sklearn.linear_model import LinearRegression

class Planeta:
    
    def __init__(self, e, a, t):
        
        self.t = t
        self.dt = t[1] - t[0] # Paso del tiempo
        
        self.e = e  # Excentricidad
        self.a_ = a # Semi-eje mayor
        
        self.G = 4*np.pi**2 # Unidades gaussianas
        
        self.r = np.zeros(3)
        self.v = np.zeros_like(self.r)
        self.a = np.zeros_like(self.r)
        
        self.r[0] = self.a_ * (1-self.e) #Periapsis distancia es mínima
        self.v[1] = np.sqrt( self.G*(1+self.e)/(self.a_*(1.-self.e)) ) #Velocidad a partir de 2 ley de Kepler
        
        self.R = np.zeros((len(t),len(self.r)))
        self.V = np.zeros_like(self.R)
        
        # El valor del pasado
        self.rp = self.r
        
    
        
    def GetAceleration(self):
        
        d = np.linalg.norm(self.r)
        self.a = -self.G/d**3*self.r
          
    def Evolution(self,i):
        
        self.SetPosition(i)
        self.SetVelocity(i)
        self.GetAceleration()
        
        if i==0:
            self.r = self.rp + self.v * self.dt
        else:
            # rp pasado, r presente rf futuro
            self.rf = 2 * self.r - self.rp + self.a * self.dt**2
            self.v = (self.rf - self.rp)/(2 * self.dt) #Derivada central
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
        for i in range(1,len(Dist)-1):
            if Dist[i] < Dist[i-1] and Dist[i] < Dist[i+1]:
                timeup.append(self.t[i])
            
        return timeup
    
def GetPlanetas(t):
    
    Mercurio = Planeta(0.2056,0.307,t)
    Venus = Planeta(0.0067,0.7233,t)
    Tierra = Planeta(0.01671,1.,t)
    Marte = Planeta(0.0933,1.524,t)
    Jupiter = Planeta(0.04877,5.2044,t)
    
    return [Mercurio,Venus,Tierra,Marte,Jupiter]

def RunSimulation(t,Planetas):
    
    for it in tqdm(range(len(t)), desc='Running simulation', unit=' Steps' ):
        #print(it)
        for i in range(len(Planetas)):
            Planetas[i].Evolution(it)       
    return Planetas

def GetPeriod(Planetas):
    periods = np.zeros(len(Planetas))
    for p in range(len(Planetas)):
        perihelio = Planetas[p].GetPerihelio()
        periods[p] = perihelio[1]-perihelio[0]
    return periods

def Geta3(Planetas):
    acubo = []
    for p in range(len(Planetas)):
        acubo.append((Planetas[p].a_) ** 3)
    return acubo
    
dt = 0.001
tmax = 30
t = np.arange(0.,tmax,dt)
Planetas_or = GetPlanetas(t)
Planetas = RunSimulation(t,Planetas_or)
periodos = GetPeriod(Planetas)
acubos = Geta3(Planetas_or)
tcuadrados = periodos ** 2
coef = np.polyfit(acubos,tcuadrados,1)
poly1d_fn = np.poly1d(coef)
names = ["Mercurio","Venus","Tierra","Marte","Jupiter"]

print("Resultados ejercicio 4")
for i in range(len(names)):
    print(names[i] + ": " + str(periodos[i]))

print("Regresión lineal")
print("Pendiente (m): " + str(coef[0]))
print("Punto de corte (b): " + str(coef[1]))
print("Masa del sol en gaussianas: "+ str(coef[0])+ " y en SI corresponde a: " + str(coef[0] * 1.989E30) + "kg")

plt.figure()
plt.scatter(acubos,tcuadrados,label = 'Simulation')
plt.plot(acubos,poly1d_fn(acubos),label = 'Linear Regression', color='r')
plt.title("Tercera Ley de Kepler")
plt.xlabel("a^3 (ua)")
plt.ylabel("T^2 (año terrestre)")
plt.legend()
plt.show()