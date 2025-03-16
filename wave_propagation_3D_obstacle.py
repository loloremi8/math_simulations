import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Simulation parameters
Lx, Ly, Lz = 2, 1, 2         # Domain size
Nx, Ny, Nz = 150, 75, 150    # Grid points
dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz
c = 150                      # Wave speed
dt = 0.00005                 # Time step
Nt = 800                     # Number of time steps

# CFL condition
if c * dt / min(dx, dy, dz) > 1:
    raise ValueError("CFL condition not satisfied. Adjust dx, dy, dz, or dt.")

# Initialize wave field
u = np.zeros((Nt, Nx, Ny, Nz))
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
z = np.linspace(0, Lz, Nz)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Initial conditions: Gaussian pulse
u[0, :, :, :] = np.exp(-100 * ((X - 0.5)**2 + (Y - 0.5)**2 + (Z - 1)**2))
u[1, :, :, :] = u[0, :, :, :]  # Initial velocity = 0

# Defining a pillar in the region
pillar_start_x = int(0.7 * Nx)  # Start of the pillar in the X direction
pillar_end_x = int(0.85 * Nx)   # End of the pillar in the X direction
pillar_start_y = int(0.3 * Ny)  # Start of the pillar in the Y direction
pillar_end_y = int(0.7 * Ny)    # End of the pillar in the Y direction
pillar_start_z = 0              # Start of the pillar in the Z direction
pillar_end_z = int(0.6 * Nz)    # End of the pillar in the Z direction

# Set the wave field values to zero in the pillar region
u[:, pillar_start_x:pillar_end_x, pillar_start_y:pillar_end_y, pillar_start_z:pillar_end_z] = 0

# Finite difference method
for n in range(1, Nt - 1):
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            for k in range(1, Nz - 1):
                u[n+1, i, j, k] = (2 * u[n, i, j, k] - u[n-1, i, j, k] +
                                   (c * dt / dx)**2 * (u[n, i+1, j, k] - 2 * u[n, i, j, k] + u[n, i-1, j, k]) +
                                   (c * dt / dy)**2 * (u[n, i, j+1, k] - 2 * u[n, i, j, k] + u[n, i, j-1, k]) +
                                   (c * dt / dz)**2 * (u[n, i, j, k+1] - 2 * u[n, i, j, k] + u[n, i, j, k-1]))
    # Ensure the pillar region remains zero
    u[n+1, pillar_start_x:pillar_end_x, pillar_start_y:pillar_end_y, pillar_start_z:pillar_end_z] = 0

# Directory to save frames
output_dir = "frames_wave_propagation_3D_obstacle"
os.makedirs(output_dir, exist_ok=True)

# Animation: Visualizing slices
fig = plt.figure(figsize=(19.2, 10.8), dpi=200)
ax = fig.add_subplot(111, projection='3d')
slice_idx = Nz // 2  # Middle slice in the z-direction
ax.view_init(elev=20, azim=80) # Set initial viewing angle

def update(frame):
    ax.clear()
    plot = ax.plot_surface(X[:, :, slice_idx], Y[:, :, slice_idx], u[frame, :, :, slice_idx], cmap="viridis")
    ax.set_title(f"3D Wave Propagation with an obstacle (Slice at z=1), Frame {frame}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u(x, y, t)")
    ax.set_zlim(-2, 2)
    ax.set_box_aspect([Lx, Ly, Lz])
    
    # Visualize the static pillar
    X_pillar, Y_pillar, Z_pillar = np.meshgrid(
        x[pillar_start_x:pillar_end_x],
        y[pillar_start_y:pillar_end_y],
        z[pillar_start_z:pillar_end_z],
        indexing='ij'
    )
    
    # Plot all six faces of the cuboid
    ax.plot_surface(X_pillar[:, :, 0], Y_pillar[:, :, 0], Z_pillar[:, :, 0], color='red', alpha=0.5)  # Bottom face
    ax.plot_surface(X_pillar[:, :, -1], Y_pillar[:, :, -1], Z_pillar[:, :, -1], color='red', alpha=0.5)  # Top face
    ax.plot_surface(X_pillar[:, 0, :], Y_pillar[:, 0, :], Z_pillar[:, 0, :], color='red', alpha=0.5)  # Front face
    ax.plot_surface(X_pillar[:, -1, :], Y_pillar[:, -1, :], Z_pillar[:, -1, :], color='red', alpha=0.5)  # Back face
    ax.plot_surface(X_pillar[0, :, :], Y_pillar[0, :, :], Z_pillar[0, :, :], color='red', alpha=0.5)  # Left face
    ax.plot_surface(X_pillar[-1, :, :], Y_pillar[-1, :, :], Z_pillar[-1, :, :], color='red', alpha=0.5)  # Right face
    
    # Save the current frame as a PNG image
    plt.savefig(f"{output_dir}/frame_{frame:04d}.png", dpi=200)
    return plot,''

# Create the animation
ani = FuncAnimation(fig, update, frames=Nt, interval=50)

# Show the plot (optional)
plt.show()