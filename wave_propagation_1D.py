import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Simulation parameters
L = 1.0       # Length of the string
c = 1.0       # Wave speed
Nx = 100      # Number of spatial points
Nt = 300      # Number of time steps
dx = L / Nx   # Spatial resolution
dt = 0.005    # Time step size

# CFL condition (stability)
if c * dt / dx > 1:
    raise ValueError("The CFL condition is not satisfied. Adjust dx or dt.")

# Initialize the wave field
u = np.zeros((Nt, Nx))
x = np.linspace(0, L, Nx)

# Initial conditions: a Gaussian pulse
u[0, :] = np.exp(-100 * (x - 0.5)**2)  # Initial displacement
u[1, :] = u[0, :]                      # Initial velocity = 0

# Finite difference method
for n in range(1, Nt - 1):
    for i in range(1, Nx - 1):
        u[n+1, i] = (2 * u[n, i] - u[n-1, i] +
                     (c * dt / dx)**2 * (u[n, i+1] - 2 * u[n, i] + u[n, i-1]))

# Animation
fig, ax = plt.subplots()
line, = ax.plot(x, u[0, :], lw=2)
ax.set_ylim(-1.2, 1.2)
ax.set_title("1D Wave Propagation")
ax.set_xlabel("x")
ax.set_ylabel("u(x, t)")

def update(frame):
    line.set_ydata(u[frame, :])
    return line,

ani = animation.FuncAnimation(fig, update, frames=Nt, interval=20)
plt.show()