import numpy as np
import matplotlib.pyplot as plt

# Henon map parameters
a = 1.4
b = 0.3

# Number of iterations
n_iterations = 140

# Initial conditions
x0, y0 = 0.1, 0.1
x0_prime, y0_prime = 0.1001, 0.1  # Slightly perturbed

# Initialize arrays to store trajectories
x, y = np.zeros(n_iterations), np.zeros(n_iterations)
x_prime, y_prime = np.zeros(n_iterations), np.zeros(n_iterations)

# Set initial values
x[0], y[0] = x0, y0
x_prime[0], y_prime[0] = x0_prime, y0_prime

# Iterate the Henon map
for n in range(1, n_iterations):
    x[n] = 1 - a * x[n - 1]**2 + y[n - 1]
    y[n] = b * x[n - 1]
    
    x_prime[n] = 1 - a * x_prime[n - 1]**2 + y_prime[n - 1]
    y_prime[n] = b * x_prime[n - 1]

# Plotting
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Phase plane plot
axs[0].plot(x, y, label="Trajectory (x0=0.1)", color='blue')
axs[0].plot(x_prime, y_prime, label="Trajectory (x0'=0.1001)", color='orange')
axs[0].set_title("Henon Map Phase Plane")
axs[0].set_xlabel("x")
axs[0].set_ylabel("y")
axs[0].legend()
axs[0].grid()

# Difference over iterations
difference = np.sqrt((x - x_prime)**2 + (y - y_prime)**2)
axs[1].plot(range(n_iterations), difference, label="Difference", color='red')
axs[1].set_title("Divergence of Trajectories Over Iterations")
axs[1].set_xlabel("Iteration")
axs[1].set_ylabel("Difference (Euclidean distance)")
axs[1].grid()
axs[1].legend()

plt.tight_layout()
plt.show()