import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Henon map parameters
a = 1.4
b = 0.3

# Number of iterations
n_iterations = 1000
n_iterations_zoom = 100

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
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
sc1 = ax1.scatter(x, y, s=1, color='blue', alpha=0.6)
ax1.set_title("Henon Map - Phase Plane")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.grid(True)

# Add a rectangle for zooming
zoom_rect = Rectangle((0, 0), 0.1, 0.1, edgecolor='red', facecolor='none')
ax1.add_patch(zoom_rect)

# Initialize the zoomed plot
sc2 = ax2.scatter([], [], s=1, color='blue', alpha=0.6)
ax2.set_title("Zoomed Henon Map")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.grid(True)

# Stack to keep track of zoom history
zoom_stack = []

def on_click(event):
    if event.inaxes != ax1:
        return
    # Update the position of the rectangle
    zoom_rect.set_xy((event.xdata, event.ydata))
    fig.canvas.draw()

def on_key(event):
    if event.key == 'enter':
        # Zoom into the rectangle area
        x_min, y_min = zoom_rect.get_xy()
        width, height = zoom_rect.get_width(), zoom_rect.get_height()
        x_max, y_max = x_min + width, y_min + height

        # Save current zoom state
        zoom_stack.append((x_min, x_max, y_min, y_max))

        # Display rendering message
        ax2.set_title("Zoomed Henon Map - Rendering...")
        fig.canvas.draw()
        plt.pause(0.1)  # Pause to force GUI update

        # Recalculate points within the zoomed area
        x_zoom = np.zeros(n_iterations_zoom)
        y_zoom = np.zeros(n_iterations_zoom)
        x_zoom[0], y_zoom[0] = (x_min + x_max) / 2, (y_min + y_max) / 2  # Start from the center of the zoom area

        print(f"Zooming into area: x_min={x_min}, x_max={x_max}, y_min={y_min}, y_max={y_max}")

        for n in range(1, n_iterations_zoom):
            x_zoom[n] = 1 - a * x_zoom[n - 1]**2 + y_zoom[n - 1]
            y_zoom[n] = b * x_zoom[n - 1]
            print(f"Iteration {n}: x={x_zoom[n]}, y={y_zoom[n]}")
            if not (x_min <= x_zoom[n] <= x_max and y_min <= y_zoom[n] <= y_max):
                print(f"Point out of bounds at iteration {n}: x={x_zoom[n]}, y={y_zoom[n]}")
                x_zoom[n] = np.nan
                y_zoom[n] = np.nan

        print(f"Zoomed points calculated: {len(x_zoom)} points")
        print(f"x_zoom: {x_zoom}")
        print(f"y_zoom: {y_zoom}")

        # Filter out NaN values
        valid_points = ~np.isnan(x_zoom) & ~np.isnan(y_zoom)
        x_zoom = x_zoom[valid_points]
        y_zoom = y_zoom[valid_points]

        print(f"Filtered x_zoom: {x_zoom}")
        print(f"Filtered y_zoom: {y_zoom}")

        # Update the scatter plot
        sc2.set_offsets(np.c_[x_zoom, y_zoom])
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y_min, y_max)
        ax2.set_title("Zoomed Henon Map")
        fig.canvas.draw()

    elif event.key == 'backspace' and zoom_stack:
        # Go back one step in the zoom
        zoom_stack.pop()
        if zoom_stack:
            x_min, x_max, y_min, y_max = zoom_stack[-1]
        else:
            x_min, x_max, y_min, y_max = ax1.get_xlim()[0], ax1.get_xlim()[1], ax1.get_ylim()[0], ax1.get_ylim()[1]

        # Recalculate points within the zoomed area
        x_zoom = np.zeros(n_iterations_zoom)
        y_zoom = np.zeros(n_iterations_zoom)
        x_zoom[0], y_zoom[0] = (x_min + x_max) / 2, (y_min + y_max) / 2  # Start from the center of the zoom area

        print(f"Reverting to area: x_min={x_min}, x_max={x_max}, y_min={y_min}, y_max={y_max}")

        for n in range(1, n_iterations_zoom):
            x_zoom[n] = 1 - a * x_zoom[n - 1]**2 + y_zoom[n - 1]
            y_zoom[n] = b * x_zoom[n - 1]
            print(f"Iteration {n}: x={x_zoom[n]}, y={y_zoom[n]}")
            if not (x_min <= x_zoom[n] <= x_max and y_min <= y_zoom[n] <= y_max):
                print(f"Point out of bounds at iteration {n}: x={x_zoom[n]}, y={y_zoom[n]}")
                x_zoom[n] = np.nan
                y_zoom[n] = np.nan

        print(f"Zoomed points recalculated: {len(x_zoom)} points")
        print(f"x_zoom: {x_zoom}")
        print(f"y_zoom: {y_zoom}")

        # Filter out NaN values
        valid_points = ~np.isnan(x_zoom) & ~np.isnan(y_zoom)
        x_zoom = x_zoom[valid_points]
        y_zoom = y_zoom[valid_points]

        print(f"Filtered x_zoom: {x_zoom}")
        print(f"Filtered y_zoom: {y_zoom}")

        # Update the scatter plot
        sc2.set_offsets(np.c_[x_zoom, y_zoom])
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y_min, y_max)
        ax2.set_title("Zoomed Henon Map")
        fig.canvas.draw()

fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()