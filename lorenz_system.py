import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Lorenz system
def lorenz(t, state, sigma, beta, rho):
    x, y, z = state
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return [dxdt, dydt, dzdt]

# Parameters
sigma = 10.0
beta = 8.0 / 3.0
rho = 28.0
initial_state = [1.0, 1.0, 1.0]  # Initial conditions
t_span = (0, 50)  # Time range
t_eval = np.linspace(*t_span, 10000)  # Time points for output

# Solve the system
solution = solve_ivp(lorenz, t_span, initial_state, args=(sigma, beta, rho), t_eval=t_eval)

# Plot the results
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot(solution.y[0], solution.y[1], solution.y[2], lw=0.5)
ax.set_title("Lorenz Attractor")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()