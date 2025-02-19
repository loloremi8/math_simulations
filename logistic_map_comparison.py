import numpy as np
import matplotlib.pyplot as plt

# Parameters
r = 3.8  # Growth rate, in the chaotic regime
n_iterations = 100  # Number of iterations
x0 = 0.5  # Initial condition 1
x0_prime = 0.5001  # Initial condition 2 (slightly perturbed)

# Logistic map function
def logistic_map(x, r):
    return r * x * (1 - x)

# Arrays to store trajectories
x = np.zeros(n_iterations)
x_prime = np.zeros(n_iterations)
x[0], x_prime[0] = x0, x0_prime

# Compute trajectories
for n in range(1, n_iterations):
    x[n] = logistic_map(x[n-1], r)
    x_prime[n] = logistic_map(x_prime[n-1], r)

# Plot results
plt.figure(figsize=(10, 6))

# Plot the two trajectories
plt.plot(range(n_iterations), x, label=f"Trajectory (x0 = {x0})", color='red')
plt.plot(range(n_iterations), x_prime, label=f"Trajectory (x0' = {x0_prime})", color='blue')

# Highlight the difference
plt.fill_between(range(n_iterations), x, x_prime, color='gray', alpha=0.5, label="Difference")

# Labels and legend
plt.title("Logistic Map: Sensitive Dependence on Initial Conditions")
plt.xlabel("Iteration (n)")
plt.ylabel("x")
plt.grid()

# Legend below the graph
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)

# Show plot
plt.tight_layout()
plt.show()