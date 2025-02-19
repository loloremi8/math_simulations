import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Lorenz system parameters
sigma = 10
rho = 28
beta = 8 / 3

# Time span
t_span = (0, 60)  # Simulate from t=0 to t=60
t_eval = np.linspace(*t_span, 6000)  # Fine time steps for smooth trajectories

# Lorenz system equations
def lorenz(t, state, sigma, rho, beta):
    x, y, z = state
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return [dxdt, dydt, dzdt]

# Initial conditions
initial_conditions_1 = [1.0, 1.0, 1.0]      # (x0, y0, z0)
initial_conditions_2 = [1.0001, 1.0, 1.0]  # Slightly perturbed

# Solve the Lorenz equations for both initial conditions
sol1 = solve_ivp(lorenz, t_span, initial_conditions_1, t_eval=t_eval, args=(sigma, rho, beta))
sol2 = solve_ivp(lorenz, t_span, initial_conditions_2, t_eval=t_eval, args=(sigma, rho, beta))

# Extract solutions
x1, y1, z1 = sol1.y
x2, y2, z2 = sol2.y

# Plotting
fig = plt.figure(figsize=(14, 6))

# 3D plot of trajectories
ax = fig.add_subplot(121, projection='3d')
ax.plot(x1, y1, z1, label="Trajectory 1 (x0=1.0)", color='blue', alpha=0.5)
ax.plot(x2, y2, z2, label="Trajectory 2 (x0=1.0001)", color='orange', alpha=0.6)
ax.set_title("Lorenz System Trajectories")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.legend()

# Difference over time
ax2 = fig.add_subplot(122)
difference = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
ax2.plot(t_eval, difference, label="Trajectory Difference", color='red', alpha=0.5)
ax2.set_title("Divergence of Trajectories Over Time")
ax2.set_xlabel("Time")
ax2.set_ylabel("Difference (Euclidean distance)")
ax2.grid()
ax2.legend()

plt.tight_layout()
plt.show()