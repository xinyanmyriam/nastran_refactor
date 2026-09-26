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
    Vector u0(n), v0(n);
    u0 << 1.0, 0.0, 0.0;
    v0 << 0.0, 0.0, 0.0;

    // Precompute constants for central difference method
    // For central difference: u_{n+1} = (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1}
    // But we need to handle the first step specially using initial velocity
    
    // Compute effective stiffness matrix for time stepping
    // K_eff = M/dt^2
    Matrix K_eff = M / (dt * dt);
    
    // Compute coefficient matrix A = (2*M/dt^2 - K)
    Matrix A = 2.0 * K_eff - K;

    // Store displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Step 0: initial displacement
    displacements.push_back(u0);

    // Compute u_{-1} using central difference formula and initial velocity
    // v0 = (u1 - u_{-1}) / (2*dt) => u_{-1} = u1 - 2*dt*v0
    // But we don't have u1 yet. Instead, use the standard central difference initialization:
    // u_{-1} = u0 - dt*v0  (first-order accurate initialization)
    Vector u_minus1 = u0 - dt * v0;

    // Now we have u_{-1}, u0, u1
    Vector u_prev = u_minus1;  // u_{n-1}
    Vector u_curr = u0;        // u_n
    Vector u_next;             // u_{n+1}

    // Perform time integration for 100 steps
    for (int step = 1; step <= nsteps; ++step) {
        // Central difference formula: u_{n+1} = 2*u_n - u_{n-1} - dt^2*K*u_n
        u_next = 2.0 * u_curr - u_prev - dt * dt * K * u_curr;
        
        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_next;
        
        // Store displacement at this step
        displacements.push_back(u_curr);
    }

    // Output JSON with exactly 11 rows: steps 0, 10, 20, ..., 100
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    
    for (int i = 0; i <= 10; ++i) {
        int step = i * 10;
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