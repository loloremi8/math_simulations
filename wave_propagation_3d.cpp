#include <vector>      // For using dynamic arrays (vectors)
#include <cmath>       // For mathematical functions like exp()
#include <iostream>    // For console input/output
#include <GL/glut.h>   // For OpenGL graphics
#include <GL/gl.h>     // For OpenGL graphics
#include <chrono>      // For high-resolution clock
#include <thread>      // For sleep function
#include <algorithm>   // For std::min and std::max

// Simulation parameters - declare as global constants
const double Lx = 1.0, Ly = 1.0, Lz = 1.0;              // Domain size
const int Nx = 200, Ny = 200, Nz = 200;                 // Resolution - grid points
const double dx = Lx / Nx, dy = Ly / Ny, dz = Lz / Nz;
const double c = 530.0;                                  // Wave speed
const double dt = 0.000005;                              // Time steps
const int Nt = 2000;                                    // Number of frames

// Global variables for visualization
std::vector<std::vector<std::vector<double>>> u_current(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
std::vector<std::vector<std::vector<double>>> u_prev(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
std::vector<std::vector<std::vector<double>>> u_next(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
int currentStep = 0;
bool simulationRunning = true;

// Function declarations
void init();
void calculateNextStep();
void display();
double calculateMaxStableC();

// OpenGL visualization functions
void init() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);  // Enable transparency
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  // Set blend function
    glClearColor(1.0, 1.0, 1.0, 1.0);  // White background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1.0, 0.1, 10.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // Adjusted camera position to view wave along x-axis
    gluLookAt(0.5, 2.5, 1.0,    // Camera position: above and behind
              0.5, -1.0, 0.0,     // Look at point: center of x-axis
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

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Draw coordinate axes
    glBegin(GL_LINES);
    glColor3f(0.3, 0.3, 0.3);  // Lighter gray for axes
    // X axis
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    // Y axis
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    // Z axis
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
    glEnd();

    // Draw grid on xy-plane - darker gray, more opaque
    glColor4f(0.2f, 0.2f, 0.2f, 0.6f);  // Darker gray, more opaque
    glLineWidth(0.7f);
    
    // Draw grid lines parallel to x-axis
    const float gridStep = 0.1f;
    for (float y = 0.0f; y <= 1.0f + gridStep/2.0f; y += gridStep) {
        glBegin(GL_LINE_STRIP);
        for (float x = 0.0f; x <= 1.0f + gridStep/2.0f; x += gridStep) {
            int i = static_cast<int>(std::min(x * Nx, (float)(Nx-1)));
            int j = static_cast<int>(std::min(y * Ny, (float)(Ny-1)));
            float height = u_current[i][j][Nz/2];
            glVertex3f(std::min(x, 1.0f), std::min(y, 1.0f), height);
        }
        glEnd();
    }
    
    // Draw grid lines parallel to y-axis
    for (float x = 0.0f; x <= 1.0f + gridStep/2.0f; x += gridStep) {
        glBegin(GL_LINE_STRIP);
        for (float y = 0.0f; y <= 1.0f + gridStep/2.0f; y += gridStep) {
            int i = static_cast<int>(std::min(x * Nx, (float)(Nx-1)));
            int j = static_cast<int>(std::min(y * Ny, (float)(Ny-1)));
            float height = u_current[i][j][Nz/2];
            glVertex3f(std::min(x, 1.0f), std::min(y, 1.0f), height);
        }
        glEnd();
    }

    // Draw wave surface with height-based coloring and transparency
    for (int i = 0; i < Nx-1; i++) {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j = 0; j < Ny; j++) {
            float height1 = u_current[i][j][Nz/2];
            float height2 = u_current[i+1][j][Nz/2];
            
            float h1 = (height1 + 1.0f) / 2.0f;
            float h2 = (height2 + 1.0f) / 2.0f;
            
            // First vertex color with alpha (transparency)
            glColor4f(
                0.4f + 0.75f * h1,  // Red
                0.0f + h1,          // Green
                1.0f - 0.8f * h1,  // Blue
                0.65f                // Alpha (0.0 = fully transparent, 1.0 = solid)
            );
            glVertex3f(i*dx, j*dy, height1);
            
            // Second vertex color with alpha
            glColor4f(
                0.3f + 0.7f * h2,
                0.0f + 1.0f * h2,
                1.0f - 0.8f * h2,
                0.65              // Alpha
            );
            glVertex3f((i+1)*dx, j*dy, height2);
        }
        glEnd();
    }
    
    // Draw contour lines - solid black, thicker
    glColor4f(0.0f, 0.0f, 0.0f, 1.0f);  // Solid black
    glLineWidth(1.2f);  // Thicker lines for emphasis
    
    const float contourInterval = 0.025f;  // Distance between contour lines
    
    for (int i = 0; i < Nx-1; i++) {
        for (int j = 0; j < Ny-1; j++) {
            float h00 = u_current[i][j][Nz/2];
            float h10 = u_current[i+1][j][Nz/2];
            float h01 = u_current[i][j+1][Nz/2];
            
            // Draw contour line if height crosses a contour level
            for (float level = -1.0f; level <= 1.0f; level += contourInterval) {
                if ((h00 <= level && h10 >= level) || 
                    (h00 >= level && h10 <= level)) {
                    // Interpolate x position where contour crosses
                    float t = (level - h00) / (h10 - h00);
                    glBegin(GL_LINES);
                    glVertex3f(i*dx + t*dx, j*dy, level);
                    glVertex3f(i*dx + t*dx, (j+1)*dy, level);
                    glEnd();
                }
                
                if ((h00 <= level && h01 >= level) || 
                    (h00 >= level && h01 <= level)) {
                    // Interpolate y position where contour crosses
                    float t = (level - h00) / (h01 - h00);
                    glBegin(GL_LINES);
                    glVertex3f(i*dx, j*dy + t*dy, level);
                    glVertex3f((i+1)*dx, j*dy + t*dy, level);
                    glEnd();
                }
            }
        }
    }
    
    glutSwapBuffers();

    if (simulationRunning) {
        calculateNextStep();
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(16));
    glutPostRedisplay();
}

double calculateMaxStableC() {
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double dz = Lz / Nz;
    
    // CFL condition in 3D
    double denominator = dt * std::sqrt(1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz));
    return 1.0 / denominator;
}

int main(int argc, char** argv) {
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
                                                    (z - 0.5)*(z - 0.5)));
                u_current[i][j][k] = u_prev[i][j][k];
            }
        }
    }

    // Initialize OpenGL
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Wave Simulation");
    init();
    glutDisplayFunc(display);
    glutMainLoop();

    return 0;
}