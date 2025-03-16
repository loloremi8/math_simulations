# Math Simulations

This repository contains various mathematical simulations implemented in Python. Each simulation demonstrates different mathematical phenomena and concepts using numerical methods and visualizations. These simulations were used in my seminar work in mathematics about dynamical systems and chaos theory.

## Simulations

### 1. Chua's Circuit Comparison
- **File:** `chua's_circuit_comparison.py`
- **Description:** Simulates and compares the trajectories of Chua's circuit for two slightly different initial conditions. Chua's circuit is a simple electronic circuit that exhibits chaotic behavior.
- **Visualization:** 3D trajectories and time series plots.

### 2. Double Pendulum Comparison
- **File:** `double_pendulum_comparison.py`
- **Description:** Simulates and compares the motion of a double pendulum for two slightly different initial conditions. The double pendulum is a classic example of a chaotic system.
- **Visualization:** Trajectories and phase space plots.

### 3. Duffing Oscillator Comparison
- **File:** `duffing_oscillator_comparison.py`
- **Description:** Simulates and compares the motion of the Duffing oscillator for two slightly different initial conditions. The Duffing oscillator is a nonlinear system that can exhibit chaotic behavior.
- **Visualization:** Time evolution and phase space plots.

### 4. Henon Map
- **File:** `henon_map.py`
- **Description:** Simulates the Henon map, a discrete-time dynamical system that exhibits chaotic behavior.
- **Visualization:** Phase plane plot.

### 5. Henon Map Comparison
- **File:** `henon_map_comparison.py`
- **Description:** Simulates and compares the trajectories of the Henon map for two slightly different initial conditions.
- **Visualization:** Phase plane plot and divergence over iterations.

### 6. Henon Map Zoom
- **File:** `henon_map_zoom.py`
- **Description:** Simulates the Henon map and allows interactive zooming into specific regions of the phase plane.
- **Visualization:** Phase plane plot with zoom functionality.
- **Warning:** This is still WIP the zoom doesn't properly  work.

### 7. Kuramoto-Sivashinsky Equation Comparison
- **File:** `kuramoto-sivashinsky_equation_comparison.py`
- **Description:** Simulates the Kuramoto-Sivashinsky equation, a nonlinear partial differential equation that exhibits chaotic behavior.
- **Visualization:** 3D surface plots of the solution and divergence.

### 8. Logistic Map
- **File:** `logistic_map.py`
- **Description:** Simulates the logistic map, a simple mathematical model of population growth that can exhibit chaotic behavior.
- **Visualization:** Bifurcation diagram.

### 9. Logistic Map Comparison
- **File:** `logistic_map_comparison.py`
- **Description:** Simulates and compares the trajectories of the logistic map for two slightly different initial conditions.
- **Visualization:** Trajectories and divergence over iterations.

### 10. Lorenz System
- **File:** `lorenz_system.py`
- **Description:** Simulates the Lorenz system, a set of differential equations that model atmospheric convection and exhibit chaotic behavior.
- **Visualization:** 3D plot of the Lorenz attractor.

### 11. Lorenz System Comparison
- **File:** `lorenz_system_comparison.py`
- **Description:** Simulates and compares the trajectories of the Lorenz system for two slightly different initial conditions.
- **Visualization:** 3D trajectories and divergence over time.

### 12. Rossler Attractor Comparison
- **File:** `rossler_attractor_comparison.py`
- **Description:** Simulates and compares the trajectories of the Rössler attractor for two slightly different initial conditions.
- **Visualization:** 3D trajectories and divergence over time.

### 13. Wave Propagation 1D
- **File:** `wave_propagation_1D.py`
- **Description:** Simulates the propagation of a wave in a 1D medium using the finite difference method.
- **Visualization:** Animation of the wave propagation.

### 14. Wave Propagation 2D
- **File:** `wave_propagation_2D.py`
- **Description:** Simulates the propagation of a wave in a 2D medium using the finite difference method.
- **Visualization:** Animation of the wave propagation.

### 15. Wave Propagation 3D
- **File:** `wave_propagation_3D.py`
- **Description:** Simulates the propagation of a wave in a 3D medium using the finite difference method.
- **Visualization:** Animation of the wave propagation.

### 16. Wave Propagation 3D with Obstacle
- **File:** `wave_propagation_3D_obstacle.py`
- **Description:** Simulates the propagation of a wave in a 3D medium with an obstacle using the finite difference method.
- **Visualization:** Animation of the wave propagation with an obstacle.
- **Warning:** This took a long time to render since it needs 1.4 billion calculations. It was render not on the faster PC, but I think I hit the problem with Python itself. For more practical use it'd be better to rewrite in C++.

### 17. Wave Propagation 3D Version 2
- **File:** `wave_propagation_3D_ver_2.py`
- **Description:** Another version of the 3D wave propagation simulation with different parameters.
- **Visualization:** Animation of the wave propagation.

### 18. Double Pendulum
- **File:** `double_pendulum.py`
- **Description:** Simulates the motion of a double pendulum.
- **Visualization:** Animation of the double pendulum motion.

### 19. GIF Maker
- **File:** `gif_maker.py`
- **Description:** Creates a GIF from a sequence of images stored in a directory.
- **Usage:** Modify the `frames_dir` and `output_filename` variables to specify the input directory and output GIF filename.

## Requirements

- Python 3.x
- NumPy
- Matplotlib
- SciPy
- Pillow (for GIF creation)

You can install all libraries with the following line:

*pip install numpy matplotlib scipy pillow*

## Usage

To run a simulation, simply run the corresponding Python file in your IDE or by the Python terminal.