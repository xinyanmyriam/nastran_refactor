#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <json/json.h> // We'll implement minimal JSON output manually to avoid external deps

// Since we can't use external JSON library, we'll write our own minimal JSON formatter
// This is a simple implementation that produces valid JSON for the required structure

std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(16) << x;
    std::string s = oss.str();
    // Remove trailing zeros but keep at least one digit after decimal
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos + 1);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    return s;
}

std::string json_array(const std::vector<double>& v) {
    std::string result = "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i > 0) result += ",";
        result += format_double(v[i]);
    }
    result += "]";
    return result;
}

std::string json_timeseries(const std::vector<std::vector<double>>& data) {
    std::string result = "[";
    for (size_t i = 0; i < data.size(); ++i) {
        if (i > 0) result += ",";
        result += json_array(data[i]);
    }
    result += "]";
    return result;
}

int main() {
    // Test case parameters
    const int n = 3; // 3-DOF system
    const double dt = 0.01;
    const int nsteps = 100;
    
    // Stiffness matrix K = [2,-1,0; -1,2,-1; 0,-1,2]
    Eigen::MatrixXd K(n, n);
    K << 2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0;
    
    // Mass matrix M = I (identity)
    Eigen::MatrixXd M(n, n);
    M.setIdentity();
    
    // Damping matrix C = 0
    Eigen::MatrixXd C(n, n);
    C.setZero();
    
    // Initial conditions
    Eigen::VectorXd u0(n), v0(n);
    u0 << 1.0, 0.0, 0.0;
    v0.setZero();
    
    // Central difference method:
    // u_{n+1} = (2*M/dt^2 - K)^(-1) * [(M/dt^2) * u_{n-1} + (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1} ???]
    // Actually, the standard central difference recurrence is:
    // M*u_{n+1} = (2*M - K*dt^2)*u_n - M*u_{n-1} + P_n*dt^2
    // But here P=0, so: u_{n+1} = (2*I - K*dt^2)*u_n - u_{n-1}
    // However, the Fortran comment says:
    //   K_eff = M/dt^2
    //   u_{n+1} = K_eff^-1 * [ (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1} ]
    // Which simplifies to: u_{n+1} = u_n - dt^2*M^{-1}*K*u_n + u_n - u_{n-1} = 2*u_n - u_{n-1} - dt^2*M^{-1}*K*u_n
    // Since M=I, this is: u_{n+1} = 2*u_n - u_{n-1} - dt^2*K*u_n
    
    // So the recurrence is: u_{n+1} = (2*I - dt^2*K)*u_n - u_{n-1}
    
    // Precompute the coefficient matrix
    Eigen::MatrixXd A = 2.0 * Eigen::MatrixXd::Identity(n, n) - dt*dt * K;
    
    // Storage for displacement history (we need steps 0,10,20,...,100 -> 11 points)
    std::vector<std::vector<double>> timeseries;
    timeseries.reserve(11);
    
    // Initialize u_{-1}, u_0, u_1
    Eigen::VectorXd u_prev(n); // u_{n-1}
    Eigen::VectorXd u_curr(n); // u_n
    Eigen::VectorXd u_next(n); // u_{n+1}
    
    // u_0 = u0
    u_curr = u0;
    
    // Compute u_{-1} using initial velocity: u_{-1} = u_0 - dt*v_0
    // But the Fortran FORM1 subroutine does: U1(I) = U0(I)-DELTT*UDOT0(I)
    u_prev = u0 - dt * v0;
    
    // Store step 0
    std::vector<double> step0;
    for (int i = 0; i < n; ++i) {
        step0.push_back(u_curr(i));
    }
    timeseries.push_back(step0);
    
    // Perform time integration
    for (int step = 1; step <= nsteps; ++step) {
        // u_{n+1} = A * u_n - u_{n-1}
        u_next = A * u_curr - u_prev;
        
        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_next;
        
        // Store every 10th step (steps 0,10,20,...,100)
        if (step % 10 == 0) {
            std::vector<double> step_data;
            for (int i = 0; i < n; ++i) {
                step_data.push_back(u_curr(i));
            }
            timeseries.push_back(step_data);
        }
    }
    
    // Output as JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":" << json_timeseries(timeseries) << "}" << std::endl;
    
    return 0;
}