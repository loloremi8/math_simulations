import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters for Chua's circuit
alpha = 9.8
beta = 14.87
m0 = -1.143
m1 = -0.714

# Nonlinearity function
def f(x):
    return m1 * x + 0.5 * (m0 - m1) * (np.abs(x + 1) - np.abs(x - 1))

# Chua's circuit equations
def chua_circuit(t, y):
    x, y_, z = y
    dxdt = alpha * (y_ - x - f(x))
    dydt = x - y_ + z
    dzdt = -beta * y_
    return [dxdt, dydt, dzdt]

# Time span for the simulation
t_span = (0, 60)
t_eval = np.linspace(t_span[0], t_span[1], 6000)

# Initial conditions
y0 = [0.1, 0.2, -0.1]        # Initial condition 1
y0_prime = [0.1001, 0.2, -0.1]  # Slight perturbation

# Solve the system for both initial conditions
sol = solve_ivp(chua_circuit, t_span, y0, t_eval=t_eval)
sol_prime = solve_ivp(chua_circuit, t_span, y0_prime, t_eval=t_eval)

# Extract solutions
x, y, z = sol.y
x_prime, y_prime, z_prime = sol_prime.y

# Visualization
fig = plt.figure(figsize=(14, 6))

# 3D plot of trajectories
ax = fig.add_subplot(121, projection='3d')
ax.plot(x, y, z, label='Initial condition 1', color='blue')
ax.plot(x_prime, y_prime, z_prime, label='Initial condition 2', color='orange')
ax.set_title("Chua's Attractor - 3D Trajectories")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")
ax.legend()

# Time series plot of x(t)
plt.subplot(1, 2, 2)
plt.plot(sol.t, x, label=r"$x(t)$, Initial: $[0.1, 0.2, -0.1]$", color="blue")
plt.plot(sol_prime.t, x_prime, label=r"$x'(t)$, Initial: $[0.1001, 0.2, -0.1]$", color="orange")
plt.title("Time Series of $x(t)$")
plt.xlabel("Time")
plt.ylabel("$x(t)$")
plt.grid(True)

# Move the legend below the plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=1)

plt.tight_layout()
plt.show()