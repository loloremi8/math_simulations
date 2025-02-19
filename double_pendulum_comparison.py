import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Constants
g = 9.81  # Acceleration due to gravity, m/s^2
l1, l2 = 1.0, 1.0  # Lengths of the pendulums, m
m1, m2 = 1.0, 1.0  # Masses of the pendulums, kg

# Equations of motion for the double pendulum
def double_pendulum_derivatives(t, y, l1, l2, m1, m2, g):
    θ1, ω1, θ2, ω2 = y
    Δθ = θ2 - θ1

    # Accelerations
    denom1 = (2 * m1 + m2 - m2 * np.cos(2 * Δθ))
    denom2 = l2 * denom1 / l1

    a1 = (-g * (2 * m1 + m2) * np.sin(θ1) - m2 * g * np.sin(θ1 - 2 * θ2) -
          2 * np.sin(Δθ) * m2 * (ω2**2 * l2 + ω1**2 * l1 * np.cos(Δθ))) / (l1 * denom1)

    a2 = (2 * np.sin(Δθ) *
          (ω1**2 * l1 * (m1 + m2) + g * (m1 + m2) * np.cos(θ1) +
           ω2**2 * l2 * m2 * np.cos(Δθ))) / denom2

    return [ω1, a1, ω2, a2]

# Initial conditions
θ1_0, θ2_0 = np.pi / 4, np.pi / 6  # Initial angles
ω1_0, ω2_0 = 0, 0  # Initial angular velocities
initial_state_1 = [θ1_0, ω1_0, θ2_0, ω2_0]

# Slightly different initial conditions
δ = 1e-3
initial_state_2 = [θ1_0 + δ, ω1_0, θ2_0, ω2_0]

# Time parameters
t_span = (0, 20)
t_eval = np.linspace(*t_span, 1000)

# Solve the differential equations
sol1 = solve_ivp(double_pendulum_derivatives, t_span, initial_state_1, t_eval=t_eval,
                 args=(l1, l2, m1, m2, g))
sol2 = solve_ivp(double_pendulum_derivatives, t_span, initial_state_2, t_eval=t_eval,
                 args=(l1, l2, m1, m2, g))

# Extract solutions
θ1_1, θ2_1 = sol1.y[0], sol1.y[2]
θ1_2, θ2_2 = sol2.y[0], sol2.y[2]

# Plot trajectories
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(t_eval, θ1_1, label=r'$\theta_1$ (Initial Condition 1)')
plt.plot(t_eval, θ1_2, label=r'$\theta_1$ (Initial Condition 2)', linestyle='--')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta_1$ (rad)')
plt.legend()
plt.title('Trajectory of $\Theta_1$')

plt.subplot(1, 2, 2)
plt.plot(θ1_1, θ2_1, label='Initial Condition 1')
plt.plot(θ1_2, θ2_2, label='Initial Condition 2', linestyle='--')
plt.xlabel(r'$\theta_1$ (rad)')
plt.ylabel(r'$\theta_2$ (rad)')
plt.legend()
plt.title('Phase Space of $\Theta_1$ and $\Theta_2$')

plt.tight_layout()
plt.show()