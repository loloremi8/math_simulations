import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Simulation parameters
Lx, Ly, Lz = 1, 1, 1         # Domain size
Nx, Ny, Nz = 40, 40, 40      # Grid points
dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz
c = 50                       # Wave speed
dt = 0.00005                 # Time step
Nt = 600                     # Number of time steps

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

# Directory to save frames
output_dir = "frames"
os.makedirs(output_dir, exist_ok=True)

# Animation: Visualizing slices
fig = plt.figure(figsize=(19.2, 10.8), dpi=200)
ax = fig.add_subplot(111, projection='3d')
slice_idx = Nz // 2  # Middle slice in the z-direction
ax.view_init(elev=15, azim=45) # Set initial viewing angle

def update(frame):
    ax.clear()
    plot = ax.plot_surface(X[:, :, slice_idx], Y[:, :, slice_idx], u[frame, :, :, slice_idx], cmap="viridis")
    ax.set_title(f"3D Wave Propagation (Slice at z=0.5), Frame {frame}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u(x, y, t)")
    ax.set_zlim(-1, 1)
    
    # Save the current frame as a PNG image
    plt.savefig(f"{output_dir}/frame_{frame:04d}.png", dpi=200)
    return plot,

# Create the animation
ani = FuncAnimation(fig, update, frames=Nt, interval=50)

# Show the plot (optional)
plt.show()