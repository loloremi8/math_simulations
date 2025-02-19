import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Simulation parameters
Lx, Ly, Lz = 1.0, 1.0, 1.0  # Domain size
Nx, Ny, Nz = 30, 30, 30     # Grid points (reduce for performance)
dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz
c = 2.5                      # Wave speed (1.0 original)
dt = 0.001                   # Time step (0.005 original)
Nt = 200                     # Number of time steps (100 original)

# CFL condition
if c * dt / min(dx, dy, dz) > 1:
    raise ValueError("CFL condition not satisfied. Adjust dx, dy, dz, or dt.")

# Initialize wave field
u = np.zeros((Nt, Nx, Ny, Nz))
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
z = np.linspace(0, Lz, Nz)
X, Y, Z = np.meshgrid(x, y, z)

# Initial conditions: Gaussian pulse
u[0, :, :, :] = np.exp(-100 * ((X - 0.5)**2 + (Y - 0.5)**2 + (Z - 0.5)**2))
u[1, :, :, :] = u[0, :, :, :]  # Initial velocity = 0

# Finite difference method
for n in range(1, Nt - 1):
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            for k in range(1, Nz - 1):
                u[n+1, i, j, k] = (2 * u[n, i, j, k] - u[n-1, i, j, k] +
                                   (c * dt / dx)**2 * (u[n, i+1, j, k] - 2 * u[n, i, j, k] + u[n, i-1, j, k]) +
                                   (c * dt / dy)**2 * (u[n, i, j+1, k] - 2 * u[n, i, j, k] + u[n, i, j-1, k]) +
                                   (c * dt / dz)**2 * (u[n, i, j, k+1] - 2 * u[n, i, j, k] + u[n, i, j, k-1]))

# Animation: Visualizing slices
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
slice_idx = Nz // 2  # Middle slice in the z-direction
plot = ax.plot_surface(X[:, :, slice_idx], Y[:, :, slice_idx], u[0, :, :, slice_idx], cmap="viridis")

ax.set_title("3D Wave Propagation (Slice at z=0.5)")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x, y, t)")
ax.set_zlim(-1, 1)

def update(frame):
    ax.clear()
    plot = ax.plot_surface(X[:, :, slice_idx], Y[:, :, slice_idx], u[frame, :, :, slice_idx], cmap="viridis")
    ax.set_title(f"3D Wave Propagation (Slice at z=0.5), Frame {frame}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u(x, y, t)")
    ax.set_zlim(-1, 1)
    return plot,

ani = FuncAnimation(fig, update, frames=Nt, interval=50)
plt.show()