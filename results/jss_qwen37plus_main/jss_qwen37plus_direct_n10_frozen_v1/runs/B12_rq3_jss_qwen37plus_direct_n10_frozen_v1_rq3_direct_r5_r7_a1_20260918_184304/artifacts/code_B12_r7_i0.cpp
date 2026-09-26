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
    // But we need to handle the first step using initial velocity
    // u1 = u0 + dt*v0 + 0.5*dt^2*a0, where a0 = M^{-1}*(P0 - C*v0 - K*u0)
    // Since C=0 and P0=0 (no external load), a0 = -M^{-1}*K*u0 = -K*u0 (since M=I)
    
    // Compute u1 using initial conditions
    Vector u1 = u0 + dt * v0 + 0.5 * dt * dt * (-K * u0);

    // Store displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);
    displacements.push_back(u0);  // step 0

    // Central difference integration
    // u_{n+1} = (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1}
    // Since M = I: u_{n+1} = (2/dt^2 - K) * u_n - (1/dt^2) * u_{n-1}
    Matrix A = 2.0 / (dt * dt) * Matrix::Identity(n, n) - K;
    Matrix B = -1.0 / (dt * dt) * Matrix::Identity(n, n);

    Vector un_minus_1 = u0;
    Vector un = u1;

    // Perform time integration
    for (int step = 1; step <= nsteps; ++step) {
        Vector un_plus_1 = A * un + B * un_minus_1;
        
        // Update for next iteration
        un_minus_1 = un;
        un = un_plus_1;
        
        if (step == nsteps) {
            displacements.push_back(un);
        } else if (step % 10 == 0) {
            // We want steps 0, 10, 20, ..., 100
            // Step 0 is already added, so add step 10, 20, ..., 100
            displacements.push_back(un);
        }
    }

    // Ensure we have exactly 11 entries (0,10,20,...,100)
    // We should have: step 0, then steps 10,20,...,100 -> total 11
    if (displacements.size() != 11) {
        // Recompute to ensure correct sampling
        displacements.clear();
        displacements.push_back(u0);  // step 0
        
        un_minus_1 = u0;
        un = u1;
        
        for (int step = 1; step <= nsteps; ++step) {
            Vector un_plus_1 = A * un + B * un_minus_1;
            un_minus_1 = un;
            un = un_plus_1;
            
            if (step % 10 == 0) {
                displacements.push_back(un);
            }
        }
    }

    // Output as JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    
    for (size_t i = 0; i < displacements.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(16) << displacements[i](j);
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;

    return 0;
}