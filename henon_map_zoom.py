import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Henon map parameters
a = 1.4
b = 0.3

# Number of iterations
n_iterations = 10000
n_iterations_zoom = 50000

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

# Plot the full Henon map in both panels initially
sc1 = ax1.scatter(x, y, s=1, color='blue', alpha=0.6)
sc2 = ax2.scatter(x, y, s=1, color='blue', alpha=0.6)  # Initialize right plot with full map

ax1.set_title("Henon Map")
ax2.set_title("Zoomed View")
for ax in [ax1, ax2]:
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)
    ax.set_xlim(-1.5, 1.5)  # Set same initial limits for both plots
    ax.set_ylim(-0.4, 0.4)

# Add a square zoom box (making sure width equals height)
rect_size = 0.2
zoom_rect = Rectangle((0, 0), rect_size, rect_size, edgecolor='red', facecolor='none')
ax1.add_patch(zoom_rect)

# Initialize zoom stack with the full view
zoom_stack = []
zoom_stack.append({
    'points': np.c_[x, y],
    'xlim': ax2.get_xlim(),
    'ylim': ax2.get_ylim(),
    'size': rect_size
})

def on_click(event):
    if event.inaxes != ax1:
        return
    # Center the square on the click
    size = zoom_rect.get_width()  # Use width since it's a square
    zoom_rect.set_xy((event.xdata - size/2, event.ydata - size/2))
    fig.canvas.draw()

def on_enter(event):
    if event.key != 'enter':
        return

    # Get the zoom area coordinates
    x_min, y_min = zoom_rect.get_xy()
    size = zoom_rect.get_width()  # Use width since it's a square
    x_max = x_min + size
    y_max = y_min + size

    # Store both points and view settings
    center_x = (x_min + x_max) / 2
    center_y = (y_min + y_max) / 2
    zoom_stack.append({
        'points': sc2.get_offsets(),
        'xlim': ax2.get_xlim(),
        'ylim': ax2.get_ylim(),
        'size': size
    })

    # Generate new points in the zoomed area
    x_zoom = np.zeros(n_iterations_zoom)
    y_zoom = np.zeros(n_iterations_zoom)
    
    # Start from multiple points near the center
    points_found = 0
    for dx in np.linspace(-size/4, size/4, 10):
        for dy in np.linspace(-size/4, size/4, 10):
            x_zoom[0] = center_x + dx
            y_zoom[0] = center_y + dy
            
            for i in range(1, n_iterations_zoom):
                x_zoom[i] = 1 - a * x_zoom[i-1]**2 + y_zoom[i-1]
                y_zoom[i] = b * x_zoom[i-1]
                
                if (x_min <= x_zoom[i] <= x_max and y_min <= y_zoom[i] <= y_max):
                    points_found += 1

            if points_found > 100:
                break
        if points_found > 100:
            break

    # Filter points within the zoom window
    mask = ((x_zoom >= x_min) & (x_zoom <= x_max) & 
            (y_zoom >= y_min) & (y_zoom <= y_max))
    x_filtered = x_zoom[mask]
    y_filtered = y_zoom[mask]

    # Update the zoomed view
    sc2.set_offsets(np.c_[x_filtered, y_filtered])
    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(y_min, y_max)
    
    # Make zoom square smaller for next zoom while keeping it centered
    new_size = size/2
    zoom_rect.set_width(new_size)
    zoom_rect.set_height(new_size)
    zoom_rect.set_xy((center_x - new_size/2, center_y - new_size/2))
    
    fig.canvas.draw()

def on_backspace(event):
    if event.key != 'backspace' or not zoom_stack:
        return
    
    # Restore previous view and points
    previous = zoom_stack.pop()
    sc2.set_offsets(previous['points'])
    ax2.set_xlim(previous['xlim'])
    ax2.set_ylim(previous['ylim'])
    
    # Restore previous zoom square size
    zoom_rect.set_width(previous['size'])
    zoom_rect.set_height(previous['size'])
    
    # Center the zoom square in the current view
    center_x = np.mean(ax2.get_xlim())
    center_y = np.mean(ax2.get_ylim())
    zoom_rect.set_xy((center_x - previous['size']/2,
                      center_y - previous['size']/2))
    
    fig.canvas.draw()

# Connect event handlers
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_enter)
fig.canvas.mpl_connect('key_press_event', on_backspace)

plt.show()