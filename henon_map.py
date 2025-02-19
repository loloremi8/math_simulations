import numpy as np
import matplotlib.pyplot as plt

# Henon map parameters
a = 1.4
b = 0.3

# Number of iterations
n_iterations = 1000

# Initial condition
x0, y0 = 0.1, 0.1

# Initialize arrays to store trajectories
x = np.zeros(n_iterations)
y = np.zeros(n_iterations)

# Set initial values
x[0], y[0] = x0, y0

# Iterate the Henon map
for n in range(1, n_iterations):
    x[n] = 1 - a * x[n - 1]**2 + y[n - 1]
    y[n] = b * x[n - 1]

# Plotting the Henon map
plt.figure(figsize=(10, 6))
plt.scatter(x, y, s=1, color='blue', alpha=0.6)
plt.title("Henon Map - Phase Plane")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()