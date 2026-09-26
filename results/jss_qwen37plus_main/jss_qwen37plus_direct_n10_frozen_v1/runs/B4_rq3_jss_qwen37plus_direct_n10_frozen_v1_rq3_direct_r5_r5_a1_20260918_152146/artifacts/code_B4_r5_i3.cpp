#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Using double precision throughout
using Real = double;
using Matrix8d = Eigen::Matrix<Real, 8, 8>;
using Vector2d = Eigen::Vector2d;
using Matrix2d = Eigen::Matrix<Real, 2, 2>;

// Helper: compute stiffness matrix for a triangular membrane element (3 nodes, 2 DOF/node => 6x6)
// Based on standard linear triangular membrane (CST) formulation
// Reference: Cook et al., "Concepts and Applications of Finite Element Analysis"
Matrix2d computeTriangularStiffness(
    const Vector2d& p1, const Vector2d& p2, const Vector2d& p3,
    Real E, Real nu, Real t) {
    
    // Material matrix D for plane stress
    Real denom = 1.0 - nu * nu;
    Matrix2d D;
    D << E / denom, E * nu / denom, E * nu / denom, E / denom;
    
    // Area of triangle
    Real area = 0.5 * std::abs((p2.x() - p1.x()) * (p3.y() - p1.y()) - 
                               (p3.x() - p1.x()) * (p2.y() - p1.y()));
    
    if (area <= 0.0) {
        return Matrix2d::Zero();
    }
    
    // Shape function derivatives (B matrix components)
    // B = [ dN1/dx  0      dN2/dx  0      dN3/dx  0
    //        0     dN1/dy    0     dN2/dy    0     dN3/dy
    //      dN1/dy dN1/dx  dN2/dy dN2/dx  dN3/dy dN3/dx ]
    // But for CST, we use constant strain triangle with:
    // N_i = a_i + b_i*x + c_i*y, where b_i, c_i are constants
    
    // Compute coefficients for shape functions
    Real x1 = p1.x(), y1 = p1.y();
    Real x2 = p2.x(), y2 = p2.y();
    Real x3 = p3.x(), y3 = p3.y();
    
    Real detJ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    if (std::abs(detJ) < 1e-15) {
        return Matrix2d::Zero();
    }
    
    // b_i and c_i coefficients
    Real b1 = y2 - y3;
    Real b2 = y3 - y1;
    Real b3 = y1 - y2;
    Real c1 = x3 - x2;
    Real c2 = x1 - x3;
    Real c3 = x2 - x1;
    
    // B matrix (3x6) for strain-displacement
    // [ b1  0  b2  0  b3  0
    //    0  c1  0  c2  0  c3
    //   c1 b1  c2 b2  c3 b3 ] / (2*area)
    Real inv2A = 1.0 / (2.0 * area);
    
    // Compute B^T * D * B * t * area
    // Since we only need the 2x2 submatrix for the in-plane DOFs at node 1 (ux,uy),
    // but actually the full 6x6 is needed. However, the original Fortran builds 4x3x3 blocks.
    // For CQUAD4 via 4 triangles, each triangle contributes to 3 nodes (6 DOFs), but only
    // the contributions to the 4 corner nodes matter.
    // So we'll compute full 6x6 for triangle and then extract/add to global 8x8.
    
    // Build full 6x6 stiffness for this triangle
    Eigen::Matrix<Real, 6, 6> K_tri = Eigen::Matrix<Real, 6, 6>::Zero();
    
    // Precompute common factor
    Real factor = t * area * inv2A * inv2A;
    
    // Loop over all 6x6 entries
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Real b_i = (i == 0) ? b1 : (i == 1) ? b2 : b3;
            Real c_i = (i == 0) ? c1 : (i == 1) ? c2 : c3;
            Real b_j = (j == 0) ? b1 : (j == 1) ? b2 : b3;
            Real c_j = (j == 0) ? c1 : (j == 1) ? c2 : c3;
            
            // K(i*2,j*2) += factor * (D(0,0)*b_i*b_j + D(0,1)*b_i*c_j + D(1,0)*c_i*b_j + D(1,1)*c_i*c_j)
            K_tri(2*i, 2*j) += factor * (D(0,0)*b_i*b_j + D(0,1)*b_i*c_j + D(1,0)*c_i*b_j + D(1,1)*c_i*c_j);
            // K(i*2,j*2+1) += factor * (D(0,0)*b_i*c_j + D(0,1)*b_i*b_j + D(1,0)*c_i*c_j + D(1,1)*c_i*b_j)
            K_tri(2*i, 2*j+1) += factor * (D(0,0)*b_i*c_j + D(0,1)*b_i*b_j + D(1,0)*c_i*c_j + D(1,1)*c_i*b_j);
            // K(i*2+1,j*2) += factor * (D(0,0)*c_i*b_j + D(0,1)*c_i*c_j + D(1,0)*b_i*b_j + D(1,1)*b_i*c_j)
            K_tri(2*i+1, 2*j) += factor * (D(0,0)*c_i*b_j + D(0,1)*c_i*c_j + D(1,0)*b_i*b_j + D(1,1)*b_i*c_j);
            // K(i*2+1,j*2+1) += factor * (D(0,0)*c_i*c_j + D(0,1)*c_i*b_j + D(1,0)*b_i*c_j + D(1,1)*b_i*b_j)
            K_tri(2*i+1, 2*j+1) += factor * (D(0,0)*c_i*c_j + D(0,1)*c_i*b_j + D(1,0)*b_i*c_j + D(1,1)*b_i*b_j);
        }
    }
    
    // Return top-left 2x2 block? No — we need full 6x6 to assemble into 8x8.
    // But the problem asks for 8x8 stiffness matrix for CQUAD4 with 4 nodes × 2 DOF.
    // So we'll build full 8x8 by assembling contributions from 4 triangles.
    // Each triangle has 3 nodes → contributes to 6 DOFs. We map local node indices to global.
    return Matrix2d::Zero(); // placeholder; full K_tri computed above
}

// Assemble 8x8 stiffness matrix for CQUAD4 using 4 triangular sub-elements
Matrix8d computeCQUAD4Stiffness(
    const std::vector<Vector2d>& nodes, // 4 nodes: (0,0), (2,0), (2,1.5), (0,1.5)
    Real E, Real nu, Real t) {
    
    // Nodes in order: A(0), B(1), C(2), D(3) as per NASTRAN convention
    // The four triangles are:
    // Triangle I:  A-B-D  (nodes[0], nodes[1], nodes[3])
    // Triangle II: B-C-A  (nodes[1], nodes[2], nodes[0])
    // Triangle III: C-D-B (nodes[2], nodes[3], nodes[1])
    // Triangle IV: D-A-C (nodes[3], nodes[0], nodes[2])
    
    // Mapping from triangle index (0-based) to local node indices (0,1,2) -> global node indices
    // From Fortran M array: M = [1,2,4,2,3,1,3,4,2,4,1,3]
    // Tri I: M1=1, M2=2, M3=4 → global nodes [0,1,3] (0-indexed)
    // Tri II: M4=2, M5=3, M6=1 → [1,2,0]
    // Tri III: M7=3, M8=4, M9=2 → [2,3,1]
    // Tri IV: M10=4, M11=1, M12=3 → [3,0,2]
    
    std::vector<std::vector<int>> triangles = {
        {0, 1, 3}, // Triangle I
        {1, 2, 0}, // Triangle II
        {2, 3, 1}, // Triangle III
        {3, 0, 2}  // Triangle IV
    };
    
    Matrix8d K = Matrix8d::Zero();
    
    // For each triangle, compute its 6x6 stiffness and add to global 8x8
    for (int tri_idx = 0; tri_idx < 4; ++tri_idx) {
        const auto& tri_nodes = triangles[tri_idx];
        const Vector2d& p1 = nodes[tri_nodes[0]];
        const Vector2d& p2 = nodes[tri_nodes[1]];
        const Vector2d& p3 = nodes[tri_nodes[2]];
        
        // Compute 6x6 stiffness for this triangle
        Eigen::Matrix<Real, 6, 6> K_tri = Eigen::Matrix<Real, 6, 6>::Zero();
        
        // Material matrix D for plane stress
        Real denom = 1.0 - nu * nu;
        Eigen::Matrix<Real, 2, 2> D;
        D << E / denom, E * nu / denom,
             E * nu / denom, E / denom;
        
        // Area of triangle
        Real area = 0.5 * std::abs((p2.x() - p1.x()) * (p3.y() - p1.y()) - 
                                   (p3.x() - p1.x()) * (p2.y() - p1.y()));
        
        if (area <= 0.0) continue;
        
        // Coefficients for shape functions
        Real x1 = p1.x(), y1 = p1.y();
        Real x2 = p2.x(), y2 = p2.y();
        Real x3 = p3.x(), y3 = p3.y();
        
        Real detJ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        if (std::abs(detJ) < 1e-15) continue;
        
        Real b1 = y2 - y3;
        Real b2 = y3 - y1;
        Real b3 = y1 - y2;
        Real c1 = x3 - x2;
        Real c2 = x1 - x3;
        Real c3 = x2 - x1;
        
        Real inv2A = 1.0 / (2.0 * area);
        Real factor = t * area * inv2A * inv2A;
        
        // Build 6x6 K_tri
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Real b_i = (i == 0) ? b1 : (i == 1) ? b2 : b3;
                Real c_i = (i == 0) ? c1 : (i == 1) ? c2 : c3;
                Real b_j = (j == 0) ? b1 : (j == 1) ? b2 : b3;
                Real c_j = (j == 0) ? c1 : (j == 1) ? c2 : c3;
                
                // K(2*i, 2*j) += factor * (D00*b_i*b_j + D01*b_i*c_j + D10*c_i*b_j + D11*c_i*c_j)
                K_tri(2*i, 2*j) += factor * (D(0,0)*b_i*b_j + D(0,1)*b_i*c_j + D(1,0)*c_i*b_j + D(1,1)*c_i*c_j);
                // K(2*i, 2*j+1) += factor * (D00*b_i*c_j + D01*b_i*b_j + D10*c_i*c_j + D11*c_i*b_j)
                K_tri(2*i, 2*j+1) += factor * (D(0,0)*b_i*c_j + D(0,1)*b_i*b_j + D(1,0)*c_i*c_j + D(1,1)*c_i*b_j);
                // K(2*i+1, 2*j) += factor * (D00*c_i*b_j + D01*c_i*c_j + D10*b_i*b_j + D11*b_i*c_j)
                K_tri(2*i+1, 2*j) += factor * (D(0,0)*c_i*b_j + D(0,1)*c_i*c_j + D(1,0)*b_i*b_j + D(1,1)*b_i*c_j);
                // K(2*i+1, 2*j+1) += factor * (D00*c_i*c_j + D01*c_i*b_j + D10*b_i*c_j + D11*b_i*b_j)
                K_tri(2*i+1, 2*j+1) += factor * (D(0,0)*c_i*c_j + D(0,1)*c_i*b_j + D(1,0)*b_i*c_j + D(1,1)*b_i*b_j);
            }
        }
        
        // Assemble into global 8x8 matrix
        // Local DOF mapping: local node i has DOFs [2*i, 2*i+1] → global DOFs [2*global_node, 2*global_node+1]
        for (int i_local = 0; i_local < 3; ++i_local) {
            int i_global = tri_nodes[i_local];
            for (int j_local = 0; j_local < 3; ++j_local) {
                int j_global = tri_nodes[j_local];
                for (int di = 0; di < 2; ++di) {
                    for (int dj = 0; dj < 2; ++dj) {
                        int i_dof = 2 * i_global + di;
                        int j_dof = 2 * j_global + dj;
                        K(i_dof, j_dof) += K_tri(2*i_local + di, 2*j_local + dj);
                    }
                }
            }
        }
    }
    
    return K;
}

// Format a double in scientific notation with fixed precision
std::string formatDouble(Real x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(12) << x;
    std::string s = oss.str();
    // Remove trailing zeros after decimal
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    return s;
}

int main() {
    // Test case inputs
    std::vector<Vector2d> nodes = {
        Vector2d(0.0, 0.0),   // Node 1
        Vector2d(2.0, 0.0),   // Node 2
        Vector2d(2.0, 1.5),   // Node 3
        Vector2d(0.0, 1.5)    // Node 4
    };
    
    Real E = 200e9;   // Pa
    Real nu = 0.3;
    Real t = 0.01;    // m
    
    // Compute 8x8 stiffness matrix
    Matrix8d K = computeCQUAD4Stiffness(nodes, E, nu, t);
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 8; ++i) {
        std::cout << "[";
        for (int j = 0; j < 8; ++j) {
            std::cout << formatDouble(K(i, j));
            if (j < 7) std::cout << ",";
        }
        std::cout << "]";
        if (i < 7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}