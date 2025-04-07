import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Henon map parameters
a = 1.4
b = 0.3

# Number of iterations
n_iterations = 20000
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
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))

# Plot the full Henon map in both panels initially
sc1 = ax1.scatter(x, y, s=1, color='blue', alpha=0.6)
sc2 = ax2.scatter(x, y, s=1, color='blue', alpha=0.6)

ax1.set_title("Henon Map")
ax2.set_title("Zoomed View")
for ax in [ax1, ax2]:
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)
    ax.set_xlim(-1.5, 1.5)
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
    size = zoom_rect.get_width()  # Use width of the square
    zoom_rect.set_xy((event.xdata - size/2, event.ydata - size/2))
    fig.canvas.draw()

def generate_points_in_bounds(x_min, x_max, y_min, y_max, n_points=None):
    """Generate points for the Henon map within specified bounds."""
    if n_points is None:
        n_points = n_iterations_zoom

    x_points = []
    y_points = []
    
    # Find a valid starting point on the attractor
    x_temp = 0.1  # Start from a known point on attractor
    y_temp = 0.1
    
    # Warmup to get onto attractor
    for _ in range(1000):
        x_new = 1 - a * x_temp**2 + y_temp
        y_new = b * x_temp
        x_temp, y_temp = x_new, y_new
    
    # Generate many more points than needed to ensure we get enough in bounds
    max_iterations = n_points * 100
    
    for _ in range(max_iterations):
        x_new = 1 - a * x_temp**2 + y_temp
        y_new = b * x_temp
        
        if (x_min <= x_new <= x_max and y_min <= y_new <= y_max):
            x_points.append(x_new)
            y_points.append(y_new)
            
            # Break if we have enough points
            if len(x_points) >= n_points:
                break
                
        x_temp, y_temp = x_new, y_new

    # Convert to numpy arrays
    x_points = np.array(x_points)
    y_points = np.array(y_points)
    
    return x_points, y_points

def on_enter(event):
    if event.key != 'enter':
        return

    # Get the zoom area coordinates
    x_min, y_min = zoom_rect.get_xy()
    size = zoom_rect.get_width()
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

    # Generate new points specifically for the zoomed area
    x_filtered, y_filtered = generate_points_in_bounds(x_min, x_max, y_min, y_max)

    # Update the zoomed view
    sc2.set_offsets(np.c_[x_filtered, y_filtered])
    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(y_min, y_max)
    
    # Make zoom square smaller for next zoom
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