#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Helper function to format double in scientific notation with fixed precision
std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos + 1);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Replace +e+ with e+
    size_t epos = s.find("+e+");
    if (epos != std::string::npos) {
        s.replace(epos, 3, "e+");
    }
    return s;
}

// Main function to compute CTRIA3 stiffness matrix
Eigen::Matrix<double, 6, 6> compute_ctria3_stiffness(
    const std::vector<std::vector<double>>& nodes,
    double E, double nu, double t) {
    
    // Nodes: 3 nodes, each with (x, y, z) - we use only x,y for 2D membrane
    // Node indexing: 0=A, 1=B, 2=C
    double x1 = nodes[0][0], y1 = nodes[0][1], z1 = nodes[0][2];
    double x2 = nodes[1][0], y2 = nodes[1][1], z2 = nodes[1][2];
    double x3 = nodes[2][0], y3 = nodes[2][1], z3 = nodes[2][2];

    // Set up the E matrix (3x3) for coordinate transformation
    Eigen::Matrix<double, 3, 3> E_mat;
    // First find I-vector = RSUBB - RSUBA (non-normalized)
    double e1 = x2 - x1;
    double e3 = y2 - y1;
    double e5 = z2 - z1;

    // Now find length = X-SUB-B
    double xsubb = std::sqrt(e1*e1 + e3*e3 + e5*e5);
    if (xsubb < 1.0e-6) {
        throw std::runtime_error("Degenerate triangle: zero length edge AB");
    }

    // Normalize I-vector
    e1 /= xsubb;
    e3 /= xsubb;
    e5 /= xsubb;

    // Store I-vector in first column of E_mat
    E_mat(0,0) = e1;
    E_mat(1,0) = e3;
    E_mat(2,0) = e5;

    // Now take RSUBC - RSUBA and store temporarily in E(2),E(4),E(6)
    double e2 = x3 - x1;
    double e4 = y3 - y1;
    double e6 = z3 - z1;

    // X-SUB-C = I . (RSUBC - RSUBA)
    double xsubc = e1*e2 + e3*e4 + e5*e6;

    // Crossing I-vector to (RSUBC-RSUBA) gives K-vector (non-normalized)
    double e7 = e3*e6 - e5*e4;
    double e8 = e5*e2 - e1*e6;
    double e9 = e1*e4 - e3*e2;

    // Length of K-vector = Y-SUB-C
    double ysubc = std::sqrt(e7*e7 + e8*e8 + e9*e9);
    if (ysubc < 1.0e-6) {
        throw std::runtime_error("Degenerate triangle: zero area");
    }

    // Normalize K-vector
    e7 /= ysubc;
    e8 /= ysubc;
    e9 /= ysubc;

    // Store K-vector in third column of E_mat
    E_mat(0,2) = e7;
    E_mat(1,2) = e8;
    E_mat(2,2) = e9;

    // J vector = K cross I
    double j1 = e5*e8 - e3*e9;
    double j3 = e1*e9 - e5*e7;
    double j5 = e3*e7 - e1*e8;

    // Normalize J-vector
    double temp = std::sqrt(j1*j1 + j3*j3 + j5*j5);
    if (temp == 0.0) {
        throw std::runtime_error("Degenerate triangle: collinear points");
    }
    j1 /= temp;
    j3 /= temp;
    j5 /= temp;

    // Store J-vector in second column of E_mat
    E_mat(0,1) = j1;
    E_mat(1,1) = j3;
    E_mat(2,1) = j5;

    // Volume of element, REELMU, FLAMDA, DELTA
    double vol = xsubb * ysubc * t / 2.0;
    double reelmu = 1.0 / xsubb;
    double flambda = 1.0 / ysubc;
    double delta = xsubc / xsubb - 1.0;

    // Form the C matrix (3x6) partitioned as CSUBA, CSUBB, CSUBC (each 3x2)
    // We'll build it as a 3x6 matrix
    Eigen::Matrix<double, 3, 6> C;
    // CSUBA (columns 0,1)
    C(0,0) = -reelmu;      C(0,1) = 0.0;
    C(1,0) = 0.0;          C(1,1) = flambda * delta;
    C(2,0) = flambda * delta; C(2,1) = -reelmu;
    // CSUBB (columns 2,3)
    C(0,2) = reelmu;       C(0,3) = 0.0;
    C(1,2) = 0.0;          C(1,3) = -flambda * reelmu * xsubc;
    C(2,2) = -flambda * reelmu * xsubc; C(2,3) = reelmu;
    // CSUBC (columns 4,5)
    C(0,4) = 0.0;          C(0,5) = 0.0;
    C(1,4) = 0.0;          C(1,5) = flambda;
    C(2,4) = flambda;      C(2,5) = 0.0;

    // Material matrix G (3x3) for isotropic plane stress
    // For isotropic material: G11 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G22 = E/(1-nu^2), G33 = E/(2*(1+nu))
    double denom = 1.0 - nu * nu;
    double g11 = E / denom;
    double g12 = nu * E / denom;
    double g22 = E / denom;
    double g33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3> G;
    G << g11, g12, 0.0,
         g12, g22, 0.0,
         0.0, 0.0, g33;

    // Compute stiffness matrix: K = VOL * (E*C)' * G * (E*C)
    // First compute E*C (3x6)
    Eigen::Matrix<double, 3, 6> EC = E_mat * C;
    
    // Then compute (EC)' * G * EC (6x6)
    Eigen::Matrix<double, 6, 6> K = vol * EC.transpose() * G * EC;
    
    return K;
}

int main() {
    // Test case: Nodes (0,0), (2,0), (1,1.5) in meters
    std::vector<std::vector<double>> nodes = {
        {0.0, 0.0, 0.0},   // Node A
        {2.0, 0.0, 0.0},   // Node B  
        {1.0, 1.5, 0.0}    // Node C
    };
    
    double E = 2.1e11;    // Pa
    double nu = 0.3;      // Poisson's ratio
    double t = 0.01;      // m

    try {
        Eigen::Matrix<double, 6, 6> K = compute_ctria3_stiffness(nodes, E, nu, t);
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 6; ++i) {
            std::cout << "[";
            for (int j = 0; j < 6; ++j) {
                std::cout << format_double(K(i,j));
                if (j < 5) std::cout << ",";
            }
            std::cout << "]";
            if (i < 5) std::cout << ",";
        }
        std::cout << "]}" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}