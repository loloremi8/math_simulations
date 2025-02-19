import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D

# Parameters for the Kuramoto-Sivashinsky equation
L = 10.0  # Domain size
N = 256    # Number of grid points
dx = L / N  # Grid spacing
dt = 0.01  # Time step
T = 5.0    # Total time

# Define the spatial domain (periodic boundary conditions)
x = np.linspace(0, L, N, endpoint=False)

# Define the initial condition (sinusoidal wave)
u0 = np.sin(2 * np.pi * x / L)

# Perturbed initial condition (slightly different from u0)
epsilon = 1e-2
u0_perturbed = np.sin(2 * np.pi * x / L) + epsilon * np.cos(4 * np.pi * x / L)

# Define the Kuramoto-Sivashinsky equation (nonlinear PDE)
def ks_rhs(t, u):
    u_xx = np.roll(u, -1) - 2*u + np.roll(u, 1)  # Second derivative
    u_xxxx = np.roll(u, -2) - 4*u + 6*u - 4*np.roll(u, 2) + np.roll(u, 3)  # Fourth derivative
    nonlinear_term = u * (np.roll(u, -1) - np.roll(u, 1)) / (2 * dx)  # Nonlinear term

    return -nonlinear_term + u_xx / dx**2 + u_xxxx / dx**4

# Define the time grid for both solutions (same for both u and u')
t_eval = np.linspace(0, T, int(T/dt))  # Time evaluation grid

# Solve using the same time grid for both the original and perturbed solutions
sol = solve_ivp(ks_rhs, [0, T], u0, t_eval=t_eval, method='RK45')
sol_perturbed = solve_ivp(ks_rhs, [0, T], u0_perturbed, t_eval=t_eval, method='RK45')

# Extract the solutions
u = sol.y
u_perturbed = sol_perturbed.y

# Create the figure and axis objects for 3D plotting
fig = plt.figure(figsize=(14, 6))

# Left plot: 3D surface plot of u(x,t) and u'(x,t)
ax1 = fig.add_subplot(121, projection='3d')

# Create the 3D surface for the original solution u(x,t)
T_grid, X_grid = np.meshgrid(t_eval, x)  # Create a meshgrid for time and space
ax1.plot_surface(X_grid, T_grid, u, cmap='inferno', edgecolor='none')
ax1.set_title("Kuramoto-Sivashinsky: u(x,t) and u'(x,t)")
ax1.set_xlabel("Space (x)")
ax1.set_ylabel("Time (t)")
ax1.set_zlabel("u(x,t)")
ax1.view_init(30, 220)

# Right plot: 3D surface plot of the divergence |u(x,t) - u'(x,t)|
ax2 = fig.add_subplot(122, projection='3d')

# Calculate the divergence between u(x,t) and u'(x,t)
divergence = np.abs(u - u_perturbed)

# Create the 3D surface for the divergence
ax2.plot_surface(X_grid, T_grid, divergence, cmap='viridis', edgecolor='none')
ax2.set_title("Divergence |u(x,t) - u'(x,t)|")
ax2.set_xlabel("Space (x)")
ax2.set_ylabel("Time (t)")
ax2.set_zlabel("Divergence")
ax2.view_init(30, 220)

# Display the plots
plt.tight_layout()
plt.show()