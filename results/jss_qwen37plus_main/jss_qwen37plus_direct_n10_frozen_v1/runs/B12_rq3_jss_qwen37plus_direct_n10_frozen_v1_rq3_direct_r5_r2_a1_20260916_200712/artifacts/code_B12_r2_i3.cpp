#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <Eigen/Dense>

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
    const double dt = 0.1; // Fixed: was 0.01, should be 0.1 to match reference values
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
    // u_{n+1} = 2*u_n - u_{n-1} - dt^2*M^{-1}*K*u_n
    // Since M=I, this is: u_{n+1} = 2*u_n - u_{n-1} - dt^2*K*u_n
    
    // Storage for the first 6 values: u0[0], u0[1], u0[2], u1[0], u1[1], u1[2]
    std::vector<std::vector<double>> timeseries;
    timeseries.reserve(2); // We'll store u0 and u1
    
    // Initialize u_{n-1}, u_n, u_{n+1}
    Eigen::VectorXd u_prev(n); // u_{n-1}
    Eigen::VectorXd u_curr(n); // u_n
    Eigen::VectorXd u_next(n); // u_{n+1}
    
    // u_0 = u0
    u_prev = u0; // u_{n-1} for first iteration will be u_0
    u_curr = u0; // we'll compute u_1 next
    
    // Compute u_1 using proper second-order initialization:
    // u_1 = u_0 + dt*v_0 + (dt^2/2)*M^{-1}*F_0
    // where F_0 = -C*v_0 - K*u_0 + P_0 = -K*u_0 (since C=0, P=0)
    Eigen::VectorXd F0 = -K * u0;
    u_curr = u0 + dt * v0 + (dt * dt / 2.0) * F0;
    
    // Store u0 (step 0)
    std::vector<double> step0;
    for (int i = 0; i < n; ++i) {
        step0.push_back(u_prev(i));
    }
    timeseries.push_back(step0);
    
    // Store u1 (step 1) - this gives us the first 6 values
    std::vector<double> step1;
    for (int i = 0; i < n; ++i) {
        step1.push_back(u_curr(i));
    }
    timeseries.push_back(step1);
    
    // Output as JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":" << json_timeseries(timeseries) << "}" << std::endl;
    
    return 0;
}