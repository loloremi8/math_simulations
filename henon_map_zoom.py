import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Henon map parameters
a = 1.4
b = 0.3

# Number of iterations
n_iterations = 1000
n_iterations_zoom = 20000

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

# Initialize the zoomed plot with the same limits as the original plot
sc2 = ax2.scatter([], [], s=1, color='blue', alpha=0.6)
ax2.set_title("Zoomed Henon Map")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylim(ax1.get_ylim())
ax2.grid(True)

# Stack to keep track of zoom history
zoom_stack = []

def on_click(event):
    if event.inaxes != ax1:
        return
    # Update the position of the rectangle
    zoom_rect.set_xy((event.xdata, event.ydata))
    fig.canvas.draw()

def on_enter(event):
    if event.key != 'enter':
        return

    # Zoom into the rectangle area
    x_min, y_min = zoom_rect.get_xy()
    width, height = zoom_rect.get_width(), zoom_rect.get_height()
    x_max, y_max = x_min + width, y_min + height

    # Save current zoom state
    zoom_stack.append((x_min, x_max, y_min, y_max, sc2.get_offsets(), zoom_rect.get_width(), zoom_rect.get_height()))

    # Display rendering message
    ax2.set_title("Zoomed Henon Map - Rendering...")
    fig.canvas.draw()
    plt.pause(0.1)  # Pause to force GUI update

    print(f"Zooming into area: x_min={x_min}, x_max={x_max}, y_min={y_min}, y_max={y_max}")

    # Recalculate points within the zoomed area
    x_zoom = np.zeros(n_iterations_zoom)
    y_zoom = np.zeros(n_iterations_zoom)
    x_zoom[0], y_zoom[0] = (x_min + x_max) / 2, (y_min + y_max) / 2  # Start from the center of the zoom area

    for n in range(1, n_iterations_zoom):
        x_zoom[n] = 1 - a * x_zoom[n - 1]**2 + y_zoom[n - 1]
        y_zoom[n] = b * x_zoom[n - 1]

    # Update the scatter plot with the zoomed Henon map
    sc2.set_offsets(np.c_[x_zoom, y_zoom])
    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(y_min, y_max)
    ax2.set_title("Zoomed Henon Map")
    fig.canvas.draw()

    # Update the zoom rectangle size
    zoom_rect.set_width(width / 2)
    zoom_rect.set_height(height / 2)
    fig.canvas.draw()

def on_backspace(event):
    if event.key != 'backspace' or not zoom_stack:
        return

    # Go back one step in the zoom
    x_min, x_max, y_min, y_max, offsets, rect_width, rect_height = zoom_stack.pop()

    print(f"Reverting to area: x_min={x_min}, x_max={x_max}, y_min={y_min}, y_max={y_max}")

    # Update the scatter plot with the previous zoom state
    sc2.set_offsets(offsets)
    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(y_min, y_max)
    ax2.set_title("Zoomed Henon Map")
    fig.canvas.draw()

    # Update the zoom rectangle size
    zoom_rect.set_width(rect_width)
    zoom_rect.set_height(rect_height)
    fig.canvas.draw()

fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_enter)
fig.canvas.mpl_connect('key_press_event', on_backspace)

plt.show()