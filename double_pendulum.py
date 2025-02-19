import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Constants
g = 9.81  # Gravity
L1, L2 = 1.0, 1.0  # Lengths of the rods
m1, m2 = 1.0, 1.0  # Masses of the pendulums
dt = 0.01  # Time step
t_max = 20  # Duration of simulation
num_frames = 1500  # Number of frames for animation

# Time array
t = np.linspace(0, t_max, num_frames)

# Equations of motion
def derivatives(state, L1, L2, m1, m2, g):
    theta1, omega1, theta2, omega2 = state
    delta = theta2 - theta1
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    den2 = L2 / L1 * den1

    dtheta1_dt = omega1
    dtheta2_dt = omega2

    domega1_dt = ((m2 * g * np.sin(theta2) * np.cos(delta) -
                   m2 * np.sin(delta) * (L2 * omega2**2 + L1 * omega1**2 * np.cos(delta)) -
                   (m1 + m2) * g * np.sin(theta1)) / den1)

    domega2_dt = (((m1 + m2) * (L1 * omega1**2 * np.sin(delta) -
                    g * np.sin(theta2) +
                    g * np.sin(theta1) * np.cos(delta)) +
                    m2 * L2 * omega2**2 * np.sin(delta) * np.cos(delta)) / den2)

    return np.array([dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt])

# Simulate motion
def simulate(theta1, theta2, omega1, omega2, L1, L2, m1, m2, g, t):
    state = np.array([theta1, omega1, theta2, omega2])
    positions = [state]
    for _ in t[1:]:
        k1 = dt * derivatives(state, L1, L2, m1, m2, g)
        k2 = dt * derivatives(state + 0.5 * k1, L1, L2, m1, m2, g)
        k3 = dt * derivatives(state + 0.5 * k2, L1, L2, m1, m2, g)
        k4 = dt * derivatives(state + k3, L1, L2, m1, m2, g)
        state = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        positions.append(state)
    return np.array(positions)

# Initial conditions
theta1_calm, theta2_calm = np.pi / 4, np.pi / 6
theta1_chaotic, theta2_chaotic = np.pi / 4 + 0.001, np.pi / 6
omega1_calm, omega2_calm = 0, 0
omega1_chaotic, omega2_chaotic = 0, 0

calm_positions = simulate(theta1_calm, theta2_calm, omega1_calm, omega2_calm, L1, L2, m1, m2, g, t)
chaotic_positions = simulate(theta1_chaotic, theta2_chaotic, omega1_chaotic, omega2_chaotic, L1, L2, m1, m2, g, t)

# Extract x, y positions
def calculate_xy(positions, L1, L2):
    theta1, _, theta2, _ = positions.T
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)
    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)
    return x1, y1, x2, y2

x1_calm, y1_calm, x2_calm, y2_calm = calculate_xy(calm_positions, L1, L2)
x1_chaotic, y1_chaotic, x2_chaotic, y2_chaotic = calculate_xy(chaotic_positions, L1, L2)

# Animation
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

for a in ax:
    a.set_xlim(-2.5, 2.5)
    a.set_ylim(-2.5, 2.5)
    a.set_aspect('equal')

lines = [ax[0].plot([], [], 'o-', lw=2)[0], ax[1].plot([], [], 'o-', lw=2)[0]]
traces = [ax[0].plot([], [], 'r-', lw=1)[0], ax[1].plot([], [], 'r-', lw=1)[0]]

ax[0].set_title("Calm Pendulum")
ax[1].set_title("Chaotic Pendulum")

def init():
    for line, trace in zip(lines, traces):
        line.set_data([], [])
        trace.set_data([], [])
    return lines + traces

def update(frame):
    # Calm pendulum
    lines[0].set_data([0, x1_calm[frame], x2_calm[frame]], [0, y1_calm[frame], y2_calm[frame]])
    traces[0].set_data(x2_calm[:frame], y2_calm[:frame])

    # Chaotic pendulum
    lines[1].set_data([0, x1_chaotic[frame], x2_chaotic[frame]], [0, y1_chaotic[frame], y2_chaotic[frame]])
    traces[1].set_data(x2_chaotic[:frame], y2_chaotic[:frame])

    return lines + traces

ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=True)

# Save as GIF
output_filename = "double_pendulum_comparison.gif"
ani.save(output_filename, writer=PillowWriter(fps=30))
print(f"Animation saved as {output_filename}")