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
    double x1 = nodes[0][0], y1 = nodes[0][1];
    double x2 = nodes[1][0], y2 = nodes[1][1];
    double x3 = nodes[2][0], y3 = nodes[2][1];

    // Compute area of triangle
    double area = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    if (area < 1.0e-12) {
        throw std::runtime_error("Degenerate triangle: zero area");
    }

    // Compute coefficients for shape function derivatives
    // For linear triangle: N_i = a_i + b_i*x + c_i*y
    // where: a_i = x_j*y_k - x_k*y_j, b_i = y_j - y_k, c_i = x_k - x_j
    // for cyclic permutations (i,j,k)
    
    // Coefficients for node 1 (A)
    double a1 = x2*y3 - x3*y2;
    double b1 = y2 - y3;
    double c1 = x3 - x2;
    
    // Coefficients for node 2 (B)
    double a2 = x3*y1 - x1*y3;
    double b2 = y3 - y1;
    double c2 = x1 - x3;
    
    // Coefficients for node 3 (C)
    double a3 = x1*y2 - x2*y1;
    double b3 = y1 - y2;
    double c3 = x2 - x1;

    // The Jacobian matrix [J] = [[b1, c1], [b2, c2], [b3, c3]] but we need [J] for the element
    // Actually, the standard approach: [B] matrix contains derivatives dN_i/dx, dN_i/dy
    // For linear triangle: dN_i/dx = b_i / (2*area), dN_i/dy = c_i / (2*area)
    
    double inv_2a = 1.0 / (2.0 * area);
    
    // Build B matrix (3x6) for strain-displacement: 
    // [B] = [[dN1/dx, 0, dN2/dx, 0, dN3/dx, 0],
    //        [0, dN1/dy, 0, dN2/dy, 0, dN3/dy],
    //        [dN1/dy, dN1/dx, dN2/dy, dN2/dx, dN3/dy, dN3/dx]]
    
    Eigen::Matrix<double, 3, 6> B;
    
    // Row 0: dN_i/dx
    B(0,0) = b1 * inv_2a; B(0,1) = 0.0;
    B(0,2) = b2 * inv_2a; B(0,3) = 0.0;
    B(0,4) = b3 * inv_2a; B(0,5) = 0.0;
    
    // Row 1: dN_i/dy  
    B(1,0) = 0.0; B(1,1) = c1 * inv_2a;
    B(1,2) = 0.0; B(1,3) = c2 * inv_2a;
    B(1,4) = 0.0; B(1,5) = c3 * inv_2a;
    
    // Row 2: dN_i/dy, dN_i/dx (for shear strain)
    B(2,0) = c1 * inv_2a; B(2,1) = b1 * inv_2a;
    B(2,2) = c2 * inv_2a; B(2,3) = b2 * inv_2a;
    B(2,4) = c3 * inv_2a; B(2,5) = b3 * inv_2a;

    // Material matrix D for isotropic plane stress
    // D = E/(1-nu^2) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
    double denom = 1.0 - nu * nu;
    double fac = E / denom;
    
    Eigen::Matrix<double, 3, 3> D;
    D << fac,      fac * nu, 0.0,
         fac * nu, fac,      0.0,
         0.0,      0.0,      fac * (1.0 - nu) / 2.0;

    // Compute stiffness matrix: K = t * area * B^T * D * B
    Eigen::Matrix<double, 6, 6> K = t * area * B.transpose() * D * B;
    
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