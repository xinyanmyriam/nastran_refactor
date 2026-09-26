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

    // Precompute initial acceleration
    // a0 = M^{-1}*(P0 - C*v0 - K*u0) = -K*u0 since P0=0, C=0, M=I
    Vector a0 = -K * u0;

    // Compute u1 using second-order Taylor expansion
    Vector u1 = u0 + dt * v0 + 0.5 * dt * dt * a0;

    // Store displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Step 0: initial displacement
    displacements.push_back(u0);

    // Step 1: displacement at first time step
    displacements.push_back(u1);

    // Initialize for time integration
    Vector u_prev = u0;   // u_{n-1}
    Vector u_curr = u1;   // u_n
    Vector u_next;         // u_{n+1}

    // Perform time integration for remaining steps (from step 2 to step 100)
    // We already have steps 0 and 1, so we need 99 more steps to reach step 100
    for (int step = 2; step <= nsteps; ++step) {
        // Central difference formula: u_{n+1} = 2*u_n - u_{n-1} + dt^2*a_n
        // where a_n = -K*u_n (since M=I, C=0, P=0)
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