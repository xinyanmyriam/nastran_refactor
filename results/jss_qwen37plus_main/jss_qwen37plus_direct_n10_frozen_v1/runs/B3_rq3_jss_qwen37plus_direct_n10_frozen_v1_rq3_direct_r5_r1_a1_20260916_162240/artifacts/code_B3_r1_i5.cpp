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
    // Extract coordinates
    Real x1 = p1(0), y1 = p1(1);
    Real x2 = p2(0), y2 = p2(1);
    Real x3 = p3(0), y3 = p3(1);
    
    // Compute area using cross product
    Real area = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    
    // Check for degenerate triangle
    constexpr Real eps = 1.0e-12;
    if (area < eps) {
        throw std::runtime_error("Degenerate triangle: zero area");
    }
    
    // Compute coefficients for shape function derivatives
    // For triangle with nodes 1,2,3, the shape functions are:
    // N1 = a1 + b1*x + c1*y, etc.
    // where b1 = y2 - y3, c1 = x3 - x2, a1 = x2*y3 - x3*y2, etc.
    
    // Coefficients for shape function derivatives
    Real b1 = y2 - y3;
    Real b2 = y3 - y1;
    Real b3 = y1 - y2;
    
    Real c1 = x3 - x2;
    Real c2 = x1 - x3;
    Real c3 = x2 - x1;
    
    // Factor for B matrix: 1/(2*area)
    Real fac = 1.0 / (2.0 * area);
    
    // Build B matrix (3x6) for plane stress
    // Rows: εxx, εyy, γxy
    // Columns: u1,v1,u2,v2,u3,v3
    Eigen::Matrix<Real, 3, 6> B;
    B.setZero();
    
    // εxx = ∂u/∂x = b1*u1 + b2*u2 + b3*u3 all times fac
    B(0,0) = b1 * fac;  // ∂N1/∂x
    B(0,2) = b2 * fac;  // ∂N2/∂x  
    B(0,4) = b3 * fac;  // ∂N3/∂x
    
    // εyy = ∂v/∂y = c1*v1 + c2*v2 + c3*v3 all times fac
    B(1,1) = c1 * fac;  // ∂N1/∂y
    B(1,3) = c2 * fac;  // ∂N2/∂y
    B(1,5) = c3 * fac;  // ∂N3/∂y
    
    // γxy = ∂u/∂y + ∂v/∂x = c1*u1 + c2*u2 + c3*u3 + b1*v1 + b2*v2 + b3*v3 all times fac
    B(2,0) = c1 * fac;  // ∂N1/∂y for u
    B(2,1) = b1 * fac;  // ∂N1/∂x for v
    B(2,2) = c2 * fac;  // ∂N2/∂y for u
    B(2,3) = b2 * fac;  // ∂N2/∂x for v
    B(2,4) = c3 * fac;  // ∂N3/∂y for u
    B(2,5) = b3 * fac;  // ∂N3/∂x for v
    
    // Material matrix D for plane stress
    // D = [D11 D12 0; D12 D22 0; 0 0 D33]
    Real D11 = E / (1.0 - nu * nu);
    Real D22 = D11;
    Real D12 = nu * E / (1.0 - nu * nu);
    Real D33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<Real, 3, 3> D;
    D << D11, D12, 0.0,
         D12, D22, 0.0,
         0.0, 0.0, D33;
    
    // Compute stiffness matrix: K = t * area * B^T * D * B
    Eigen::Matrix<Real, 6, 6> K = t * area * B.transpose() * D * B;
    
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