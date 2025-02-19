import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Simulation parameters
Lx, Ly = 1.0, 1.0  # Domain size
Nx, Ny = 100, 100  # Number of points in x and y
dx, dy = Lx / Nx, Ly / Ny  # Spatial resolution
c = 1.0             # Wave speed
dt = 0.005          # Time step
Nt = 300            # Number of time steps

# CFL condition (stability)
if c * dt / min(dx, dy) > 1:
    raise ValueError("CFL condition not satisfied. Adjust dx, dy, or dt.")

# Initialize the wave field
u = np.zeros((Nt, Nx, Ny))
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y)

# Initial conditions: a Gaussian pulse
u[0, :, :] = np.exp(-100 * ((X - 0.5)**2 + (Y - 0.5)**2))
u[1, :, :] = u[0, :, :]  # Initial velocity = 0

# Finite difference method
for n in range(1, Nt - 1):
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            u[n+1, i, j] = (2 * u[n, i, j] - u[n-1, i, j] +
                            (c * dt / dx)**2 * (u[n, i+1, j] - 2 * u[n, i, j] + u[n, i-1, j]) +
                            (c * dt / dy)**2 * (u[n, i, j+1] - 2 * u[n, i, j] + u[n, i, j-1]))

# Animation
fig, ax = plt.subplots()
cax = ax.imshow(u[0, :, :], extent=[0, Lx, 0, Ly], cmap="viridis", vmin=-1, vmax=1)
fig.colorbar(cax)
ax.set_title("2D Wave Propagation")
ax.set_xlabel("x")
ax.set_ylabel("y")

def update(frame):
    cax.set_array(u[frame, :, :])
    return cax,

ani = FuncAnimation(fig, update, frames=Nt, interval=20)
plt.show()