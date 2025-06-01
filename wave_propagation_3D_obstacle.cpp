#include <vector>
#include <cmath>
#include <iostream>
#include <GL/glut.h>
#include <GL/gl.h>
#include <chrono>
#include <thread>
#include <algorithm>

// Simulation parameters
const double Lx = 2.0, Ly = 1.0, Lz = 2.0;              // Domain size
const int Nx = 160, Ny = 80, Nz = 160;                  // Resolution - grid points
const double dx = Lx / Nx, dy = Ly / Ny, dz = Lz / Nz;
const double c = 1000.0;                                  // Wave speed
const double dt = 0.000005;                              // Time steps
const int Nt = 2000;                                     // Number of frames

// Pillar parameters
const int pillar_start_x = static_cast<int>(0.7 * Nx);
const int pillar_end_x = static_cast<int>(0.85 * Nx);
const int pillar_start_y = static_cast<int>(0.3 * Ny);
const int pillar_end_y = static_cast<int>(0.7 * Ny);
const int pillar_start_z = 0;
const int pillar_end_z = static_cast<int>(0.6 * Nz);

// Global variables
std::vector<std::vector<std::vector<double>>> u_current(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
std::vector<std::vector<std::vector<double>>> u_prev(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
std::vector<std::vector<std::vector<double>>> u_next(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
int currentStep = 0;
bool simulationRunning = true;

double calculateMaxStableC() {
    // CFL condition in 3D
    double denominator = dt * std::sqrt(1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz));
    return 1.0 / denominator;
}

void init() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 2.0, 0.1, 10.0);  // Modified aspect ratio for wider view
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // Adjusted camera position for better view of the elongated domain
    gluLookAt(1.0, 2.5, 1.25,    // Camera position: above and behind
              1.0, -1.0, 0.0,    // Look at point: center of x-axis
              0.0, 0.0, 1.0);    // Up vector: aligned with z-axis
}

void calculateNextStep() {
    if (currentStep >= Nt) {
        simulationRunning = false;
        return;
    }

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                // Skip calculation if point is inside pillar
                if (i >= pillar_start_x && i < pillar_end_x &&
                    j >= pillar_start_y && j < pillar_end_y &&
                    k >= pillar_start_z && k < pillar_end_z) {
                    u_next[i][j][k] = 0;
                    continue;
                }

                u_next[i][j][k] = 2.0 * u_current[i][j][k] - u_prev[i][j][k] +
                    (c * dt / dx) * (c * dt / dx) * (u_current[i+1][j][k] - 2.0 * u_current[i][j][k] + u_current[i-1][j][k]) +
                    (c * dt / dy) * (c * dt / dy) * (u_current[i][j+1][k] - 2.0 * u_current[i][j][k] + u_current[i][j-1][k]) +
                    (c * dt / dz) * (c * dt / dz) * (u_current[i][j][k+1] - 2.0 * u_current[i][j][k] + u_current[i][j][k-1]);
            }
        }
    }

    u_prev = u_current;
    u_current = u_next;
    currentStep++;
}

void drawPillar() {
    glColor4f(0.8f, 0.2f, 0.2f, 0.5f);  // Semi-transparent red
    
    float x1 = pillar_start_x * dx;
    float x2 = pillar_end_x * dx;
    float y1 = pillar_start_y * dy;
    float y2 = pillar_end_y * dy;
    float z1 = pillar_start_z * dz;
    float z2 = pillar_end_z * dz;

    // Draw the pillar as a solid cuboid
    glBegin(GL_QUADS);
    // Front face
    glVertex3f(x1, y1, z1);
    glVertex3f(x2, y1, z1);
    glVertex3f(x2, y1, z2);
    glVertex3f(x1, y1, z2);
    
    // Back face
    glVertex3f(x1, y2, z1);
    glVertex3f(x2, y2, z1);
    glVertex3f(x2, y2, z2);
    glVertex3f(x1, y2, z2);
    
    // Right face
    glVertex3f(x2, y1, z1);
    glVertex3f(x2, y2, z1);
    glVertex3f(x2, y2, z2);
    glVertex3f(x2, y1, z2);
    
    // Left face
    glVertex3f(x1, y1, z1);
    glVertex3f(x1, y2, z1);
    glVertex3f(x1, y2, z2);
    glVertex3f(x1, y1, z2);
    
    // Top face
    glVertex3f(x1, y1, z2);
    glVertex3f(x2, y1, z2);
    glVertex3f(x2, y2, z2);
    glVertex3f(x1, y2, z2);
    
    // Bottom face
    glVertex3f(x1, y1, z1);
    glVertex3f(x2, y1, z1);
    glVertex3f(x2, y2, z1);
    glVertex3f(x1, y2, z1);
    glEnd();
}

// Add these helper functions after the global variables:
struct Color {
    float r, g, b;
};

Color getHeightColor(float height) {
    // Normalize height to [-1,1] range and apply non-linear scaling
    float h = height;
    
    // Apply non-linear scaling to make middle range more sensitive
    if (std::abs(h) <= 0.4f) {
        // Scale middle range (-0.4 to 0.4) to full range
        h = h * 1.25f;
    } else {
        // Compress outer ranges
        h = (h > 0) ? 0.5f + (h - 0.4f) * 0.5f : -0.5f + (h + 0.4f) * 0.5f;
    }
    
    // Now normalize to [0,1] range
    h = (h + 1.0f) / 2.0f;
    
    // Define color gradient
    if (h < 0.2f) {
        // Deep blue to blue
        return {0.0f, 0.0f, 0.3f + 3.5f * h};
    } else if (h < 0.4f) {
        // Blue to cyan
        float t = (h - 0.2f) * 5.0f;
        return {0.0f, t, 1.0f};
    } else if (h < 0.6f) {
        // Cyan to green
        float t = (h - 0.4f) * 5.0f;
        return {0.0f, 1.0f, 1.0f - t};
    } else if (h < 0.8f) {
        // Green to yellow
        float t = (h - 0.6f) * 5.0f;
        return {t, 1.0f, 0.0f};
    } else {
        // Yellow to red
        float t = (h - 0.8f) * 5.0f;
        return {1.0f, 1.0f - t, 0.0f};
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Draw coordinate axes
    glBegin(GL_LINES);
    glColor3f(0.3, 0.3, 0.3);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(2.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 2.0);
    glEnd();

    // Draw grid lines
    glColor4f(0.2f, 0.2f, 0.2f, 0.6f);
    glLineWidth(0.7f);
    
    // Draw grid lines parallel to x-axis
    const float gridStep = 0.1f;
    for (float y = 0.0f; y <= 1.0f + gridStep/2.0f; y += gridStep) {
        glBegin(GL_LINE_STRIP);
        for (float x = 0.0f; x <= 2.0f + gridStep/2.0f; x += gridStep) {
            int i = static_cast<int>(std::min(x * Nx / 2, (float)(Nx-1)));
            int j = static_cast<int>(std::min(y * Ny, (float)(Ny-1)));
            float height = u_current[i][j][Nz/2];
            glVertex3f(std::min(x, 2.0f), std::min(y, 1.0f), height);
        }
        glEnd();
    }
    
    // Draw grid lines parallel to y-axis
    for (float x = 0.0f; x <= 2.0f + gridStep/2.0f; x += gridStep) {
        glBegin(GL_LINE_STRIP);
        for (float y = 0.0f; y <= 1.0f + gridStep/2.0f; y += gridStep) {
            int i = static_cast<int>(std::min(x * Nx / 2, (float)(Nx-1)));
            int j = static_cast<int>(std::min(y * Ny, (float)(Ny-1)));
            float height = u_current[i][j][Nz/2];
            glVertex3f(std::min(x, 2.0f), std::min(y, 1.0f), height);
        }
        glEnd();
    }

    // Draw wave surface with enhanced color gradient
    glLineWidth(1.0f);
    for (int i = 0; i < Nx-1; i++) {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j = 0; j < Ny; j++) {
            float height1 = u_current[i][j][Nz/2];
            float height2 = u_current[i+1][j][Nz/2];
            
            // Get colors for both vertices
            Color c1 = getHeightColor(height1);
            Color c2 = getHeightColor(height2);
            
            // First vertex
            glColor4f(c1.r, c1.g, c1.b, 0.75f);
            glVertex3f(i*dx, j*dy, height1);
            
            // Second vertex
            glColor4f(c2.r, c2.g, c2.b, 0.75f);
            glVertex3f((i+1)*dx, j*dy, height2);
        }
        glEnd();
    }

    // Draw contour lines with darker coloring
    glLineWidth(1.5f);  // Made lines slightly thicker
    
    const float contourInterval = 0.05f;
    for (int i = 0; i < Nx-1; i++) {
        for (int j = 0; j < Ny-1; j++) {
            float h00 = u_current[i][j][Nz/2];
            float h10 = u_current[i+1][j][Nz/2];
            float h01 = u_current[i][j+1][Nz/2];
            
            for (float level = -1.0f; level <= 1.0f; level += contourInterval) {
                // Get base color and darken it
                Color levelColor = getHeightColor(level);
                float darkFactor = 0.5f;  // Make contour lines darker
                glColor4f(
                    levelColor.r * darkFactor, 
                    levelColor.g * darkFactor, 
                    levelColor.b * darkFactor, 
                    0.9f  // More opaque
                );
                
                if ((h00 <= level && h10 >= level) || 
                    (h00 >= level && h10 <= level)) {
                    float t = (level - h00) / (h10 - h00);
                    glBegin(GL_LINES);
                    glVertex3f(i*dx + t*dx, j*dy, level);
                    glVertex3f(i*dx + t*dx, (j+1)*dy, level);
                    glEnd();
                }
                
                if ((h00 <= level && h01 >= level) || 
                    (h00 >= level && h01 <= level)) {
                    float t = (level - h00) / (h01 - h00);
                    glBegin(GL_LINES);
                    glVertex3f(i*dx, j*dy + t*dy, level);
                    glVertex3f((i+1)*dx, j*dy + t*dy, level);
                    glEnd();
                }
            }
        }
    }

    // Draw the pillar
    drawPillar();
    
    glutSwapBuffers();

    if (simulationRunning) {
        calculateNextStep();
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(16));
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    // Check CFL condition
    double maxC = calculateMaxStableC();
    std::cout << "Maximum stable wave speed (c): " << maxC << std::endl;
    if (c > maxC) {
        std::cerr << "Warning: Current c (" << c << ") exceeds maximum stable value!" << std::endl;
        std::cerr << "Simulation may become unstable." << std::endl;
    }

    // Initialize wave field with Gaussian pulse
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;
                u_prev[i][j][k] = std::exp(-100.0 * (
                    (x - 0.5)*(x - 0.5) +
                    (y - 0.5)*(y - 0.5) + 
                    (z - 1.0)*(z - 1.0)));
                u_current[i][j][k] = u_prev[i][j][k];

                // Set initial conditions to zero inside the pillar
                if (i >= pillar_start_x && i < pillar_end_x &&
                    j >= pillar_start_y && j < pillar_end_y &&
                    k >= pillar_start_z && k < pillar_end_z) {
                    u_prev[i][j][k] = 0;
                    u_current[i][j][k] = 0;
                }
            }
        }
    }

    // Initialize OpenGL
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1200, 600);  // Wider window for elongated domain
    glutCreateWindow("Wave Simulation with Obstacle");
    init();
    glutDisplayFunc(display);
    glutMainLoop();

    return 0;
}