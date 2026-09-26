#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Main rod element stiffness computation
Eigen::Matrix<double, 12, 12> compute_crod_stiffness() {
    // Test case parameters
    const double xA = 0.0, yA = 0.0, zA = 0.0;
    const double xB = 2.0, yB = 0.0, zB = 0.0;
    const double E = 200.0e9;      // Pa
    const double A = 0.01;         // m^2
    const double G = 76.923e9;     // Pa
    const double J = 5.0e-6;       // m^4
    
    // Compute length
    const double dx = xB - xA;
    const double dy = yB - yA;
    const double dz = zB - zA;
    const double L = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // Direction cosines
    const double nx = dx / L;
    const double ny = dy / L;
    const double nz = dz / L;
    
    // Build N matrix (3x3): outer product of direction vector
    Eigen::Matrix<double, 3, 3> N;
    N << nx*nx, nx*ny, nx*nz,
         nx*ny, ny*ny, ny*nz,
         nx*nz, ny*nz, nz*nz;
    
    // Compute axial and torsional stiffness coefficients
    const double dscl = E * A / L;  // axial stiffness coefficient
    const double dscr = G * J / L;  // torsional stiffness coefficient
    
    // Build full 12x12 stiffness matrix
    Eigen::Matrix<double, 12, 12> K_global = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Axial stiffness contribution (affects ux, uy, uz DOFs at both nodes: indices 0-2 and 6-8)
    K_global.block<3,3>(0,0)   =  dscl * N;
    K_global.block<3,3>(0,6)   = -dscl * N;
    K_global.block<3,3>(6,0)   = -dscl * N;
    K_global.block<3,3>(6,6)   =  dscl * N;
    
    // Torsional stiffness contribution (affects rx, ry, rz DOFs at both nodes: indices 3-5 and 9-11)
    K_global.block<3,3>(3,3)   =  dscr * N;
    K_global.block<3,3>(3,9)   = -dscr * N;
    K_global.block<3,3>(9,3)   = -dscr * N;
    K_global.block<3,3>(9,9)   =  dscr * N;
    
    return K_global;
}

int main() {
    // Compute the 12x12 stiffness matrix
    Eigen::Matrix<double, 12, 12> K = compute_crod_stiffness();
    
    // Extract top-left 6x6 block: K(A,A) — DOFs [ux,uy,uz,rx,ry,rz] of node A
    Eigen::Matrix<double, 6, 6> K_6x6 = K.block<6,6>(0,0);
    
    // Output as JSON to stdout
    std::cout << "{\n";
    std::cout << "  \"stiffness_matrix\": [\n";

    for (int i = 0; i < 6; ++i) {
        std::cout << "    [";
        for (int j = 0; j < 6; ++j) {
            std::cout << std::scientific << std::setprecision(16) << K_6x6(i, j);
            if (j < 5) std::cout << ", ";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
        std::cout << "\n";
    }

    std::cout << "  ]\n";
    std::cout << "}\n";

    return 0;
}