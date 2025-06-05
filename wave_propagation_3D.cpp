#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <GL/glut.h>   // GLUT must be included before gl.h
#include <GL/gl.h>
#include <GL/freeglut.h>  // Add this for better text rendering
#include <chrono>      // For high-resolution clock
#include <thread>      // For sleep function
#include <algorithm>   // For std::min and std::max
#include <filesystem>  // For file system operations
#include <fstream>    // For file output

// Simulation parameters - declare as global constants
const double Lx = 1.0, Ly = 1.0, Lz = 1.0;              // Domain size
const int Nx = 180, Ny = 180, Nz = 180;                 // Resolution - grid points
const double dx = Lx / Nx, dy = Ly / Ny, dz = Lz / Nz;
const double c = 310.0;                                  // Wave speed
const double dt = 0.00001;                              // Time steps
const int Nt = 1000;                                    // Number of frames

// Global variables for visualization
std::vector<std::vector<std::vector<double>>> u_current(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
std::vector<std::vector<std::vector<double>>> u_prev(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
std::vector<std::vector<std::vector<double>>> u_next(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
int currentStep = 0;
bool simulationRunning = true;

// Add after other global variables
int frameCount = 0;
std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;
float fps = 0.0f;

// In the global variables section:
//const std::string outputDir = "wave_frames";
//bool recording = true;  // Control recording

// Function declarations
void init();
void calculateNextStep();
void display();
double calculateMaxStableC();
//void saveFrame(int frameNumber);  // Commented out function declaration

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

    #pragma omp parallel for collapse(3)
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
    // Calculate FPS
    auto currentTime = std::chrono::high_resolution_clock::now();
    float deltaTime = std::chrono::duration<float>(currentTime - lastTime).count();
    if (deltaTime >= 1.0f) {
        fps = frameCount / deltaTime;
        frameCount = 0;
        lastTime = currentTime;
        // Removed console output
    }
    frameCount++;

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
    glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
    glLineWidth(1.2f);
    
    const float contourInterval = 0.025f;
    for (int i = 0; i < Nx-1; i++) {
        for (int j = 0; j < Ny-1; j++) {
            float h00 = u_current[i][j][Nz/2];
            float h10 = u_current[i+1][j][Nz/2];
            float h01 = u_current[i][j+1][Nz/2];
            
            // Draw contour lines
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
    
    // Draw frame counter and FPS
    glDisable(GL_DEPTH_TEST);  // Disable depth testing for 2D text
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 800, 0, 800, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // Background box - made smaller to match HELVETICA_10 font size
    glColor4f(0.0f, 0.0f, 0.0f, 0.7f);
    glBegin(GL_QUADS);
    glVertex2f(5, 770);    // Moved up
    glVertex2f(200, 770);  // Made narrower (from 250 to 200)
    glVertex2f(200, 790);  // Top position unchanged
    glVertex2f(5, 790);
    glEnd();

    // Text rendering - slightly adjusted position
    glColor3f(1.0f, 1.0f, 1.0f);
    glRasterPos2f(10.0f, 778.0f);  // Adjusted y position (from 775 to 778)
    
    std::stringstream ss;
    ss << "Frame: " << currentStep << "/" << Nt 
       << " | FPS: " << static_cast<int>(fps)
       << " | " << static_cast<int>(100.0f * currentStep / Nt) << "%";
    
    std::string text = ss.str();
    for (const char& c : text) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, c);
    }

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glEnable(GL_DEPTH_TEST);  // Re-enable depth testing

    if (simulationRunning) {
        calculateNextStep();
    }

    // Remove fixed sleep delay and use variable timing
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

//void saveFrame(int frameNumber) {
//    // Create output directory if it doesn't exist
//    std::filesystem::create_directories(outputDir);
//    
//    // Prepare filename
//    std::stringstream ss;
//    ss << outputDir << "/frame_" << std::setw(5) << std::setfill('0') << frameNumber << ".ppm";
//    std::string filename = ss.str();
//    
//    // Get the window size
//    GLint viewport[4];
//    glGetIntegerv(GL_VIEWPORT, viewport);
//    int width = viewport[2];
//    int height = viewport[3];
//    
//    // Allocate memory for the pixel data
//    std::vector<unsigned char> pixels(3 * width * height);
//    
//    // Read pixels from framebuffer
//    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
//    
//    // Save as PPM file
//    std::ofstream out(filename, std::ios::binary);
//    out << "P6\n" << width << " " << height << "\n255\n";
//    
//    // Flip the image vertically while writing
//    for (int y = height - 1; y >= 0; y--) {
//        out.write(reinterpret_cast<char*>(pixels.data() + y * width * 3), width * 3);
//    }
//}

int main(int argc, char** argv) {
    // Initialize timer for FPS calculation
    lastTime = std::chrono::high_resolution_clock::now();
    
    // Check wave speed stability
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