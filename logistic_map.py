import numpy as np
import matplotlib.pyplot as plt

# Parameters
r_min, r_max = 2.5, 4.0   # Range of r values to explore
n_iterations = 1000       # Number of iterations for each r
n_transient = 200         # Number of transient iterations to discard
x0 = 0.5                  # Initial condition for x

# Range of r values
r_values = np.linspace(r_min, r_max, 1000)
x_values = np.zeros((len(r_values), n_iterations - n_transient))

# Iterate for each r
for i, r in enumerate(r_values):
    x = x0
    for n in range(n_iterations):
        x = r * x * (1 - x)
        if n >= n_transient:  # Record values after transient phase
            x_values[i, n - n_transient] = x

# Plotting the bifurcation diagram
plt.figure(figsize=(10, 6))
for i, r in enumerate(r_values):
    plt.plot([r] * (n_iterations - n_transient), x_values[i, :], ',k', alpha=0.25)

plt.title("Bifurcation Diagram of the Logistic Map")
plt.xlabel("Growth Rate (r)")
plt.ylabel("Population (x)")
plt.grid(alpha=0.4)
plt.show()