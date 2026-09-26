#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Use double precision throughout
using Real = double;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

int main() {
    // Test case parameters
    const int n = 3;  // 3-DOF system
    const Real dt = 0.01;
    const int nsteps = 100;

    // Stiffness matrix K = [2,-1,0; -1,2,-1; 0,-1,2]
    Matrix K(n, n);
    K << 2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0;

    // Mass matrix M = I (identity)
    Matrix M = Matrix::Identity(n, n);

    // Damping matrix C = 0
    Matrix C = Matrix::Zero(n, n);

    // Initial conditions
    Vector u0(n);
    u0 << 1.0, 0.0, 0.0;
    
    Vector v0(n);
    v0 << 0.0, 0.0, 0.0;

    // Precompute constants for central difference method
    // K_eff = M/dt^2
    Matrix K_eff = M / (dt * dt);
    
    // For central difference: u_{n+1} = K_eff^-1 * [(2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1}]
    // So: u_{n+1} = K_eff^-1 * [A * u_n - B * u_{n-1}] where A = 2*M/dt^2 - K, B = M/dt^2
    Matrix A = 2.0 * M / (dt * dt) - K;
    Matrix B = M / (dt * dt);

    // Solve for K_eff inverse (since K_eff is diagonal and positive definite, we can use direct inverse)
    // But for generality and to match NASTRAN's approach, we'll compute the inverse
    Matrix K_eff_inv = K_eff.inverse();

    // Store displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Step 0: initial displacement
    displacements.push_back(u0);

    // Compute u_{-1} using initial velocity: u_{-1} = u_0 - dt * v_0
    Vector u_minus1 = u0 - dt * v0;

    // Current and previous displacements
    Vector u_prev = u_minus1;  // u_{n-1}
    Vector u_curr = u0;        // u_n

    // Perform time integration using central difference
    for (int step = 1; step <= nsteps; ++step) {
        // u_{n+1} = K_eff^-1 * [A * u_n - B * u_{n-1}]
        Vector u_next = K_eff_inv * (A * u_curr - B * u_prev);
        
        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_next;
        
        displacements.push_back(u_curr);
    }

    // Output JSON with exactly 11 rows: steps 0, 10, 20, ..., 100
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    
    for (int i = 0; i <= 10; ++i) {
        int step = i * 10;  // 0, 10, 20, ..., 100
        const Vector& u = displacements[step];
        
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(16) << u(j);
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;

    return 0;
}