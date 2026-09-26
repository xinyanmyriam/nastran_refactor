#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Using double precision throughout
using Real = double;

// Helper function to format double in scientific notation with fixed precision
std::string to_sci(Real x, int precision = 6) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(precision) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Replace +e+ with e+
    size_t epos = s.find("+e+");
    if (epos != std::string::npos) {
        s.replace(epos, 3, "e+");
    }
    epos = s.find("+e-");
    if (epos != std::string::npos) {
        s.replace(epos, 3, "e-");
    }
    return s;
}

// JSON-safe string escaping (minimal for numbers)
std::string json_escape(const std::string& s) {
    return s;
}

// Main computation function for CTRIA3 stiffness matrix
Eigen::Matrix<Real, 6, 6> compute_ctria3_stiffness(
    const Eigen::Vector2d& p1,  // node 1 (x,y)
    const Eigen::Vector2d& p2,  // node 2 (x,y)
    const Eigen::Vector2d& p3,  // node 3 (x,y)
    Real E,                     // Young's modulus
    Real nu,                    // Poisson's ratio
    Real t                       // thickness
) {
    // Constants
    constexpr Real eps = 1.0e-6;
    
    // Extract coordinates (z=0 for all)
    Real x1 = p1(0), y1 = p1(1), z1 = 0.0;
    Real x2 = p2(0), y2 = p2(1), z2 = 0.0;
    Real x3 = p3(0), y3 = p3(1), z3 = 0.0;
    
    // Compute E matrix (3x3): columns are i, j, k vectors of local coordinate system
    Eigen::Matrix<Real, 3, 3> E_mat;
    
    // I-vector = p2 - p1
    Real e1 = x2 - x1;
    Real e3 = y2 - y1;
    Real e5 = z2 - z1;
    
    Real xsubb = std::sqrt(e1*e1 + e3*e3 + e5*e5);
    if (xsubb < eps) {
        throw std::runtime_error("Degenerate triangle: zero length edge between nodes 1 and 2");
    }
    
    // Normalize I-vector
    e1 /= xsubb;
    e3 /= xsubb;
    e5 /= xsubb;
    
    // Store temporarily: E(2),E(4),E(6) = p3 - p1
    Real e2 = x3 - x1;
    Real e4 = y3 - y1;
    Real e6 = z3 - z1;
    
    // XSUBC = I . (p3-p1)
    Real xsubc = e1*e2 + e3*e4 + e5*e6;
    
    // K-vector = I × (p3-p1) (non-normalized)
    Real e7 = e3*e6 - e5*e4;
    Real e8 = e5*e2 - e1*e6;
    Real e9 = e1*e4 - e3*e2;
    
    Real ysubc = std::sqrt(e7*e7 + e8*e8 + e9*e9);
    if (ysubc < eps) {
        throw std::runtime_error("Degenerate triangle: collinear nodes");
    }
    
    // Normalize K-vector
    e7 /= ysubc;
    e8 /= ysubc;
    e9 /= ysubc;
    
    // J-vector = K × I
    Real j1 = e5*e8 - e3*e9;
    Real j3 = e1*e9 - e5*e7;
    Real j5 = e3*e7 - e1*e8;
    
    Real temp = std::sqrt(j1*j1 + j3*j3 + j5*j5);
    if (temp < eps) {
        throw std::runtime_error("Degenerate triangle: zero J-vector");
    }
    
    j1 /= temp;
    j3 /= temp;
    j5 /= temp;
    
    // Fill E matrix: columns are i, j, k
    E_mat << e1, j1, e7,
             e3, j3, e8,
             e5, j5, e9;
    
    // Volume of element
    Real vol = xsubb * ysubc * t / 2.0;
    
    // Material constants
    Real reelmu = 1.0 / xsubb;
    Real flambda = 1.0 / ysubc;
    Real delta = xsubc / xsubb - 1.0;
    
    // Build C matrix (3x6) partitioned as [C_A | C_B | C_C], each 3x2
    // C is stored column-wise in Fortran order, but we'll build as row-major for clarity
    Eigen::Matrix<Real, 3, 6> C;
    C.setZero();
    
    // C_A (columns 0,1)
    C(0,0) = -reelmu;   C(0,1) = 0.0;
    C(1,0) = 0.0;       C(1,1) = flambda*delta;
    C(2,0) = 0.0;       C(2,1) = -reelmu;
    
    // C_B (columns 2,3)
    C(0,2) = reelmu;    C(0,3) = 0.0;
    C(1,2) = 0.0;       C(1,3) = -flambda*reelmu*xsubc;
    C(2,2) = 0.0;       C(2,3) = reelmu;
    
    // C_C (columns 4,5)
    C(0,4) = 0.0;       C(0,5) = 0.0;
    C(1,4) = 0.0;       C(1,5) = flambda;
    C(2,4) = 0.0;       C(2,5) = flambda;
    
    // Material matrix G (3x3) for isotropic plane stress
    // G = [G11 G12 0; G12 G22 0; 0 0 G33]
    Real G11 = E / (1.0 - nu*nu);
    Real G22 = G11;
    Real G12 = nu * E / (1.0 - nu*nu);
    Real G33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<Real, 3, 3> G;
    G << G11, G12, 0.0,
         G12, G22, 0.0,
         0.0, 0.0, G33;
    
    // Compute stiffness: K = vol * (E * C)^T * G * (E * C)
    // First compute E*C: (3x3) * (3x6) = (3x6)
    Eigen::Matrix<Real, 3, 6> EC = E_mat * C;
    
    // Then compute (EC)^T * G * EC: (6x3) * (3x3) * (3x6) = (6x6)
    Eigen::Matrix<Real, 6, 6> K = vol * EC.transpose() * G * EC;
    
    return K;
}

int main() {
    try {
        // Test case: Nodes at (0,0), (2,0), (1,1.5) meters
        Eigen::Vector2d p1(0.0, 0.0);
        Eigen::Vector2d p2(2.0, 0.0);
        Eigen::Vector2d p3(1.0, 1.5);
        
        // Material properties
        Real E = 2.1e11;   // Pa
        Real nu = 0.3;      // Poisson's ratio
        Real t = 0.01;      // m
        
        // Compute stiffness matrix
        Eigen::Matrix<Real, 6, 6> K = compute_ctria3_stiffness(p1, p2, p3, E, nu, t);
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 6; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << "[";
            for (int j = 0; j < 6; ++j) {
                if (j > 0) std::cout << ",";
                std::cout << to_sci(K(i,j));
            }
            std::cout << "]";
        }
        std::cout << "]}" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}