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

    // Precompute dt squared
    const Real dt2 = dt * dt;

    // Store displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Step 0: initial displacement
    displacements.push_back(u0);

    // Compute u_{-1} using initial velocity and acceleration: 
    // u_{-1} = u_0 - dt * v_0 + 0.5 * dt^2 * a_0
    // where a_0 = M^{-1} * (-K * u_0) = -K * u_0 (since M = I)
    Vector a0 = -K * u0;
    Vector u_minus1 = u0 - dt * v0 + 0.5 * dt2 * a0;

    // Current and previous displacements
    Vector u_prev = u_minus1;  // u_{n-1}
    Vector u_curr = u0;        // u_n

    // Perform time integration using central difference
    // u_{n+1} = (2I - dt^2*K)*u_n - u_{n-1}  (since M = I, C = 0)
    for (int step = 1; step <= nsteps; ++step) {
        // u_{n+1} = (2I - dt^2*K)*u_n - u_{n-1}
        Vector u_next = 2.0 * u_curr - dt2 * K * u_curr - u_prev;
        
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