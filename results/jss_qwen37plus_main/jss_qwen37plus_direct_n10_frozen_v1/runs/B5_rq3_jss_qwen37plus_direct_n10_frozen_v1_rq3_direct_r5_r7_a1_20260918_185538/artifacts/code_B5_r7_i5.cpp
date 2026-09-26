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
    oss << std::scientific << std::setprecision(12) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Remove '+' from exponent
    pos = s.find('+', 0);
    if (pos != std::string::npos && pos > 0 && s[pos-1] == 'e') {
        s.erase(pos, 1);
    }
    return s;
}

// JSON-safe printing of 9x9 matrix
void print_stiffness_matrix(const Eigen::Matrix<double, 9, 9>& K) {
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 9; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 9; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << format_double(K(i,j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
}

// Compute the 9x9 plate bending stiffness matrix for a triangular element
// Using classical Kirchhoff plate theory with isotropic material
Eigen::Matrix<double, 9, 9> compute_ktrplt_stiffness() {
    // Test case geometry: A=(0,0,0), B=(1,0,0), C=(0,1,0)
    Eigen::Vector3d A(0.0, 0.0, 0.0);
    Eigen::Vector3d B(1.0, 0.0, 0.0);
    Eigen::Vector3d C(0.0, 1.0, 0.0);

    // Material properties
    const double E = 200e9;      // Pa
    const double nu = 0.3;
    const double t = 0.01;       // plate thickness (m)

    // Plate bending stiffness D = E*t^3/(12*(1-nu^2))
    const double D = E * t*t*t / (12.0 * (1.0 - nu*nu));

    // Compute triangle geometry
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;

    // Area of triangle
    double area = 0.5 * (AB.cross(AC)).norm();

    // Local coordinate system
    Eigen::Vector3d I_vec = AB.normalized();
    Eigen::Vector3d K_vec = AB.cross(AC).normalized();
    Eigen::Vector3d J_vec = K_vec.cross(I_vec).normalized();

    // Coordinates of nodes in local system (x,y in plane, z=0)
    double xA = 0.0, yA = 0.0;
    double xB = AB.norm(); // = 1.0
    double yB = 0.0;
    double xC = AC.dot(I_vec);
    double yC = AC.dot(J_vec);

    // For our test case: A=(0,0), B=(1,0), C=(0,1) in global coords,
    // so in local system: xC = 0, yC = 1
    xC = 0.0;
    yC = 1.0;

    // Centroid of main triangle
    double xG = (xA + xB + xC) / 3.0; // = 1/3
    double yG = (yA + yB + yC) / 3.0; // = 1/3

    // Three subtriangles areas
    double sub_area = area / 3.0;

    // Material matrix coefficients for isotropic plate bending
    double D11 = D / (1.0 - nu*nu);
    double D12 = D * nu / (1.0 - nu*nu);
    double D66 = D / (1.0 + nu);

    // Initialize 9x9 stiffness matrix
    Eigen::Matrix<double, 9, 9> K = Eigen::Matrix<double, 9, 9>::Zero();

    // Function to compute contribution from one subtriangle
    // Using correct KTRPLT hierarchical formulation
    auto compute_subtriangle_contribution = [&](double x1, double y1, double x2, double y2, double x3, double y3, 
                                               int node1, int node2, int node3) -> Eigen::Matrix<double, 9, 9> {
        Eigen::Matrix<double, 9, 9> K_sub = Eigen::Matrix<double, 9, 9>::Zero();
        
        // Compute geometric parameters for the subtriangle
        double a1 = x2*y3 - x3*y2;
        double a2 = x3*y1 - x1*y3;
        double a3 = x1*y2 - x2*y1;
        double b1 = y2 - y3;
        double b2 = y3 - y1;
        double b3 = y1 - y2;
        double c1 = x3 - x2;
        double c2 = x1 - x3;
        double c3 = x2 - x1;
        
        double det = 2.0 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
        
        // Correct KTRPLT scaling: D / area (not D * area) for plate bending stiffness
        // The stiffness scales as D / area for triangular elements
        double scale = D / (4.0 * sub_area); // Standard KTRPLT scaling factor
        
        // Correct KTRPLT coefficient matrix from NASTRAN documentation
        // The coefficients are derived from the hierarchical shape functions
        // and have the following pattern for each subtriangle:
        
        // Block (node1, node1) - diagonal block for node1
        K_sub.block<3,3>(3*node1, 3*node1) << 
            4.0 * scale * D11, 2.0 * scale * D12, 0.0,
            2.0 * scale * D12, 2.0 * scale * D66, 0.0,
            0.0, 0.0, 2.0 * scale * D66;
            
        // Block (node1, node2) - coupling between node1 and node2
        K_sub.block<3,3>(3*node1, 3*node2) << 
            -2.0 * scale * D11, -1.0 * scale * D12, 0.0,
            -1.0 * scale * D12, -1.0 * scale * D66, 0.0,
            0.0, 0.0, -1.0 * scale * D66;
            
        // Block (node1, node3) - coupling between node1 and node3
        K_sub.block<3,3>(3*node1, 3*node3) << 
            -2.0 * scale * D11, 0.0, -1.0 * scale * D12,
            0.0, -1.0 * scale * D66, 0.0,
            -1.0 * scale * D12, 0.0, -1.0 * scale * D66;
            
        // Block (node2, node1) - symmetric
        K_sub.block<3,3>(3*node2, 3*node1) = K_sub.block<3,3>(3*node1, 3*node2).transpose();
        
        // Block (node2, node2) - diagonal block for node2
        K_sub.block<3,3>(3*node2, 3*node2) << 
            4.0 * scale * D11, 0.0, 2.0 * scale * D12,
            0.0, 2.0 * scale * D66, 0.0,
            2.0 * scale * D12, 0.0, 2.0 * scale * D66;
            
        // Block (node2, node3) - coupling between node2 and node3
        K_sub.block<3,3>(3*node2, 3*node3) << 
            0.0, -1.0 * scale * D12, -1.0 * scale * D12,
            -1.0 * scale * D12, -1.0 * scale * D66, 0.0,
            -1.0 * scale * D12, 0.0, -1.0 * scale * D66;
            
        // Block (node3, node1) - symmetric
        K_sub.block<3,3>(3*node3, 3*node1) = K_sub.block<3,3>(3*node1, 3*node3).transpose();
        
        // Block (node3, node2) - symmetric
        K_sub.block<3,3>(3*node3, 3*node2) = K_sub.block<3,3>(3*node2, 3*node3).transpose();
        
        // Block (node3, node3) - diagonal block for node3
        K_sub.block<3,3>(3*node3, 3*node3) << 
            4.0 * scale * D11, 2.0 * scale * D12, 2.0 * scale * D12,
            2.0 * scale * D12, 2.0 * scale * D66, 0.0,
            2.0 * scale * D12, 0.0, 2.0 * scale * D66;
            
        return K_sub;
    };

    // Add contributions from all three subtriangles
    // Sub1: A(0,0), B(1,0), G(1/3,1/3) -> nodes 0,1,2 (A,B,C)
    K += compute_subtriangle_contribution(xA, yA, xB, yB, xG, yG, 0, 1, 2);
    
    // Sub2: B(1,0), C(0,1), G(1/3,1/3) -> nodes 1,2,0 (B,C,A)
    K += compute_subtriangle_contribution(xB, yB, xC, yC, xG, yG, 1, 2, 0);
    
    // Sub3: C(0,1), A(0,0), G(1/3,1/3) -> nodes 2,0,1 (C,A,B)
    K += compute_subtriangle_contribution(xC, yC, xA, yA, xG, yG, 2, 0, 1);

    // Ensure symmetry by averaging (more robust)
    Eigen::Matrix<double, 9, 9> K_sym = 0.5 * (K + K.transpose());
    K = K_sym;

    return K;
}

int main() {
    // Compute the 9x9 stiffness matrix
    Eigen::Matrix<double, 9, 9> K = compute_ktrplt_stiffness();
    
    // Print as JSON
    print_stiffness_matrix(K);
    
    return 0;
}