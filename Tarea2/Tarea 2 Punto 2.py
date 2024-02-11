import numpy as np
import matplotlib.pyplot as plt

def verlet(x0, v0, w, num_steps, h):
    x = np.zeros(num_steps)
    v = np.zeros(num_steps)
    x[0] = x0
    v[0] = v0

    for i in range(1, num_steps):
        x[i] = x[i-1] + v[i-1] * h + 0.5 * (-w**2) * x[i-1] * h**2
        v[i] = v[i-1] + 0.5 * (-w**2) * (x[i-1] + x[i]) * h

    return x, v

x0 = 1.0
v0 = 1.0
w = np.pi
num_steps = 1000
h = 0.1

x,v = verlet(x0, v0, w, num_steps, h)

t = np.arange(0, num_steps * h, h)
plt.plot(t, x, label='Position')
plt.plot(t, v, label='Velocity')
plt.title('Verlet Algorithm - Harmonic Oscillator Simulation')
plt.xlabel('Time')
plt.legend()
plt.show()