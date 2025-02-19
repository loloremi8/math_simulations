import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters for the Rössler attractor
a, b, c = 0.2, 0.2, 5.7

# Define the Rössler equations
def rossler(t, state):
    x, y, z = state
    dxdt = -y - z
    dydt = x + a * y
    dzdt = b + z * (x - c)
    return [dxdt, dydt, dzdt]

# Initial conditions
initial_conditions_1 = [1, 1, 1]
initial_conditions_2 = [1.001, 1, 1]

# Time range
t_span = (0, 100)
t_eval = np.linspace(t_span[0], t_span[1], 10000)

# Solve the differential equations
solution_1 = solve_ivp(rossler, t_span, initial_conditions_1, t_eval=t_eval, method='RK45')
solution_2 = solve_ivp(rossler, t_span, initial_conditions_2, t_eval=t_eval, method='RK45')

# Extract solutions
x1, y1, z1 = solution_1.y
x2, y2, z2 = solution_2.y

# Calculate divergence
divergence = np.abs(x1 - x2)

# Create a single figure with two subplots
fig = plt.figure(figsize=(16, 6))

# 3D Plot of the Rössler attractor with both trajectories (SDIC visualization)
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot(x1, y1, z1, color='blue', lw=0.5, label="Trajectory: (1, 1, 1)")
ax1.plot(x2, y2, z2, color='red', lw=0.5, linestyle='dashed', label="Trajectory: (1.001, 1, 1)")
ax1.set_title("Rössler Attractor with SDIC Visualization")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
ax1.legend()

# 2D Plot of the divergence comparison
ax2 = fig.add_subplot(122)
ax2.plot(t_eval, divergence, label="|x(t) - x'(t)|", color='green')
ax2.set_title("Divergence Between Trajectories")
ax2.set_xlabel("Time (t)")
ax2.set_ylabel("Divergence")
ax2.legend()
ax2.grid()

# Adjust layout and show the combined plots
plt.tight_layout()
plt.show()