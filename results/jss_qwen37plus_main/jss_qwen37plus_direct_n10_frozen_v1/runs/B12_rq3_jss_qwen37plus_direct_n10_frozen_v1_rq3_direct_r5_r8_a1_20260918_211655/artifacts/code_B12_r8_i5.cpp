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
    // For central difference: A * u_{n+1} = 2*K_eff*u_n - K_eff*u_{n-1}
    // where K_eff = M/dt^2 and A = 2*K_eff - K
    Matrix K_eff = M / (dt * dt);
    Matrix A = 2.0 * K_eff - K;

    // Store displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Step 0: initial displacement
    displacements.push_back(u0);

    // Compute u_{-1} using central difference formula and initial velocity
    // u_{-1} = u0 - dt*v0 + 0.5*dt^2*a0, where a0 = M^{-1}*(P0 - C*v0 - K*u0)
    // With P0=0, C=0, M=I: a0 = -K*u0
    Vector a0 = -K * u0;
    Vector u_minus1 = u0 - dt * v0 + 0.5 * dt * dt * a0;

    // Now we have u_{-1}, u0, u1
    Vector u_prev = u_minus1;  // u_{n-1}
    Vector u_curr = u0;        // u_n
    Vector u_next;             // u_{n+1}

    // Perform time integration for 100 steps
    for (int step = 1; step <= nsteps; ++step) {
        // Central difference formula: A * u_{n+1} = 2*K_eff*u_n - K_eff*u_{n-1}
        // Solve for u_{n+1} using LDLT decomposition for numerical stability
        Vector rhs = 2.0 * K_eff * u_curr - K_eff * u_prev;
        u_next = A.ldlt().solve(rhs);
        
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