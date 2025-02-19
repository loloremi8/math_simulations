import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Duffing oscillator parameters
delta = 0.2  # Damping coefficient
alpha = -1   # Linear stiffness coefficient
beta = 1     # Nonlinear stiffness coefficient
gamma = 0.3  # Forcing amplitude
omega = 1    # Forcing frequency

# Time range for simulation
t_span = (0, 50)  # Start and end time
t_eval = np.linspace(t_span[0], t_span[1], 1000)  # Time points for output

# Duffing oscillator differential equation
def duffing_oscillator(t, y, delta, alpha, beta, gamma, omega):
    x, v = y  # Position (x) and velocity (v)
    dxdt = v
    dvdt = -delta * v - alpha * x - beta * x**3 + gamma * np.cos(omega * t)
    return [dxdt, dvdt]

# Initial conditions
y0 = [0.1, 0]      # x(0) = 0.1, x'(0) = 0
y0_prime = [0.1001, 0]  # Slightly perturbed: x(0) = 0.1001, x'(0) = 0

# Solve the system for both initial conditions
sol = solve_ivp(duffing_oscillator, t_span, y0, args=(delta, alpha, beta, gamma, omega), t_eval=t_eval)
sol_prime = solve_ivp(duffing_oscillator, t_span, y0_prime, args=(delta, alpha, beta, gamma, omega), t_eval=t_eval)

# Extract solutions
x, v = sol.y
x_prime, v_prime = sol_prime.y

# Plot the results
plt.figure(figsize=(14, 6))

# Plot x(t) for both initial conditions
plt.subplot(1, 2, 1)
plt.plot(sol.t, x, label=r"$x(t)$, initial $x(0) = 0.1$", color="blue")
plt.plot(sol_prime.t, x_prime, label=r"$x'(t)$, initial $x(0) = 0.1001$", color="orange")
plt.title("Time Evolution of $x(t)$")
plt.xlabel("Time")
plt.ylabel("$x(t)$")
plt.grid(True)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)

# Phase space plot (x vs. dx/dt)
plt.subplot(1, 2, 2)
plt.plot(x, v, label=r"$x$ vs $dx/dt$ for $x(0) = 0.1$", color="blue")
plt.plot(x_prime, v_prime, label=r"$x'$ vs $dx'/dt$ for $x'(0) = 0.1001$", color="orange")
plt.title("Phase Space")
plt.xlabel("$x$")
plt.ylabel("$dx/dt$")
plt.grid(True)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)

plt.tight_layout()
plt.show()