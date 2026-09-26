#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Using double precision throughout
using Real = double;
using Matrix8x8 = Eigen::Matrix<Real, 8, 8>;
using Vector2 = Eigen::Vector2d;
using Matrix2x2 = Eigen::Matrix2d;

// Helper function to compute triangular membrane stiffness (KTRMEM equivalent)
// For isotropic linear elastic membrane: 2D plane stress
// Returns 3x3 local stiffness for triangle with nodes A,B,C (each 2-DOF: ux,uy)
Eigen::Matrix<Real, 9, 1> ktrmem_tri(const Vector2& A, const Vector2& B_point, const Vector2& C,
                                      Real E, Real nu, Real t, Real cos_th, Real sin_th) {
    // Material matrix for plane stress
    Real D11 = E / (1.0 - nu*nu);
    Real D12 = nu * D11;
    Real D22 = D11;
    Real D66 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<Real, 3, 3> D;
    D << D11, D12, 0.0,
         D12, D22, 0.0,
         0.0, 0.0, D66;
    
    // Triangle area
    Real area = std::abs((B_point.x() - A.x())*(C.y() - A.y()) - (C.x() - A.x())*(B_point.y() - A.y())) * 0.5;
    if (area == 0.0) {
        Eigen::Matrix<Real, 9, 1> K;
        K.setZero();
        return K;
    }
    
    // Strain-displacement matrix B (3x6)
    // B = [ dN1/dx  0     dN2/dx  0     dN3/dx  0
    //        0    dN1/dy   0    dN2/dy   0    dN3/dy
    //      dN1/dy dN1/dx dN2/dy dN2/dx dN3/dy dN3/dx ]
    
    // Shape function derivatives (constant for linear triangle)
    Real detJ = 2.0 * area; // Jacobian determinant
    
    // dN1/dx = (y2-y3)/detJ, dN1/dy = (x3-x2)/detJ, etc.
    Real y2_y3 = B_point.y() - C.y();
    Real x3_x2 = C.x() - B_point.x();
    Real y3_y1 = C.y() - A.y();
    Real x1_x3 = A.x() - C.x();
    Real y1_y2 = A.y() - B_point.y();
    Real x2_x1 = B_point.x() - A.x();
    
    Real dN1dx = y2_y3 / detJ;
    Real dN1dy = x3_x2 / detJ;
    Real dN2dx = y3_y1 / detJ;
    Real dN2dy = x1_x3 / detJ;
    Real dN3dx = y1_y2 / detJ;
    Real dN3dy = x2_x1 / detJ;
    
    Eigen::Matrix<Real, 3, 6> B_matrix;
    B_matrix << dN1dx, 0.0,     dN2dx, 0.0,     dN3dx, 0.0,
                0.0,     dN1dy, 0.0,     dN2dy, 0.0,     dN3dy,
                dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx;
    
    // Rotate B matrix if needed (for anisotropic or rotated material)
    // For isotropic, rotation doesn't change D, but we still apply coordinate transformation
    // Here we assume the material orientation is given by cos_th, sin_th
    // We'll build a 2x2 rotation matrix and extend to 3x3 for strain transformation
    Eigen::Matrix<Real, 2, 2> R;
    R << cos_th, -sin_th,
         sin_th,  cos_th;
    
    // Strain transformation matrix T (3x3) for plane stress
    Eigen::Matrix<Real, 3, 3> T;
    T << R(0,0)*R(0,0), R(0,1)*R(0,1), 2.0*R(0,0)*R(0,1),
         R(1,0)*R(1,0), R(1,1)*R(1,1), 2.0*R(1,0)*R(1,1),
         R(0,0)*R(1,0), R(0,1)*R(1,1), R(0,0)*R(1,1) + R(0,1)*R(1,0);
    
    // Transformed material matrix
    Eigen::Matrix<Real, 3, 3> D_rot = T * D * T.transpose();
    
    // Element stiffness in local coordinates: K = t * B^T * D * B * area
    // But since we have constant B, it's: K = t * area * B^T * D * B
    Eigen::Matrix<Real, 6, 6> K_local = t * area * B_matrix.transpose() * D_rot * B_matrix;
    
    // Return full 6x6 as vector for now
    Eigen::Matrix<Real, 36, 1> K_full;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K_full(i*6 + j) = K_local(i, j);
        }
    }
    
    return Eigen::Matrix<Real, 9, 1>::Zero(); // placeholder
}

// Standard CQUAD4 stiffness matrix computation (isoparametric, 2x2 Gauss)
Matrix8x8 compute_cquad4_stiffness(const std::vector<Vector2>& nodes,
                                   Real E, Real nu, Real t) {
    // Nodes: 4 corners in order: (0,0), (2,0), (2,1.5), (0,1.5)
    // So node 0: (0,0), node 1: (2,0), node 2: (2,1.5), node 3: (0,1.5)
    
    // Material matrix for plane stress
    Real D11 = E / (1.0 - nu*nu);
    Real D12 = nu * D11;
    Real D22 = D11;
    Real D66 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<Real, 3, 3> D;
    D << D11, D12, 0.0,
         D12, D22, 0.0,
         0.0, 0.0, D66;
    
    // Gauss points for 2x2 integration
    std::vector<Real> xi = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    std::vector<Real> eta = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    Real w = 1.0; // weight for each point in 2x2 Gauss
    
    Matrix8x8 K = Matrix8x8::Zero();
    
    // Shape functions for bilinear quadrilateral
    auto N = [](Real xi, Real eta) -> std::array<Real, 4> {
        return {
            0.25 * (1.0 - xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 + eta),
            0.25 * (1.0 - xi) * (1.0 + eta)
        };
    };
    
    // Derivatives of shape functions w.r.t. xi, eta
    auto dN_dxi = [](Real xi, Real eta) -> std::array<Real, 4> {
        return {
            -0.25 * (1.0 - eta),
            0.25 * (1.0 - eta),
            0.25 * (1.0 + eta),
            -0.25 * (1.0 + eta)
        };
    };
    
    auto dN_deta = [](Real xi, Real eta) -> std::array<Real, 4> {
        return {
            -0.25 * (1.0 - xi),
            -0.25 * (1.0 + xi),
            0.25 * (1.0 + xi),
            0.25 * (1.0 - xi)
        };
    };
    
    // Loop over Gauss points
    for (Real xi_val : xi) {
        for (Real eta_val : eta) {
            // Compute shape function derivatives w.r.t. x,y via Jacobian
            auto dN_xi = dN_dxi(xi_val, eta_val);
            auto dN_eta = dN_deta(xi_val, eta_val);
            
            // Compute Jacobian J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
            Eigen::Matrix<Real, 2, 2> J = Eigen::Matrix<Real, 2, 2>::Zero();
            for (int i = 0; i < 4; ++i) {
                J(0,0) += dN_xi[i] * nodes[i].x();
                J(0,1) += dN_xi[i] * nodes[i].y();
                J(1,0) += dN_eta[i] * nodes[i].x();
                J(1,1) += dN_eta[i] * nodes[i].y();
            }
            
            // Inverse Jacobian
            Real detJ = J.determinant();
            if (std::abs(detJ) < 1e-15) continue;
            
            Eigen::Matrix<Real, 2, 2> Jinv = J.inverse();
            
            // Compute derivatives w.r.t. x,y
            std::array<Eigen::Vector2d, 4> dN_dx;
            for (int i = 0; i < 4; ++i) {
                dN_dx[i] << Jinv(0,0) * dN_xi[i] + Jinv(0,1) * dN_eta[i],
                            Jinv(1,0) * dN_xi[i] + Jinv(1,1) * dN_eta[i];
            }
            
            // Strain-displacement matrix B (3x8)
            // B = [ dN1/dx  0     dN2/dx  0     dN3/dx  0     dN4/dx  0
            //        0    dN1/dy   0    dN2/dy   0    dN3/dy   0    dN4/dy
            //      dN1/dy dN1/dx dN2/dy dN2/dx dN3/dy dN3/dx dN4/dy dN4/dx ]
            Eigen::Matrix<Real, 3, 8> B = Eigen::Matrix<Real, 3, 8>::Zero();
            for (int i = 0; i < 4; ++i) {
                int col = i * 2;
                B(0, col) = dN_dx[i].x();     // dNi/dx
                B(1, col+1) = dN_dx[i].y();   // dNi/dy
                B(2, col) = dN_dx[i].y();     // dNi/dy
                B(2, col+1) = dN_dx[i].x();   // dNi/dx
            }
            
            // Integrand: B^T * D * B * detJ * w
            Eigen::Matrix<Real, 8, 8> K_integrand = B.transpose() * D * B * detJ * w;
            K += K_integrand;
        }
    }
    
    // Multiply by thickness
    K *= t;
    
    return K;
}

int main() {
    // Test case data
    std::vector<Vector2> nodes = {
        Vector2(0.0, 0.0),    // node 0
        Vector2(2.0, 0.0),    // node 1
        Vector2(2.0, 1.5),    // node 2
        Vector2(0.0, 1.5)     // node 3
    };
    
    Real E = 200e9;   // Pa
    Real nu = 0.3;
    Real t = 0.01;    // m
    
    // Compute stiffness matrix
    Matrix8x8 K = compute_cquad4_stiffness(nodes, E, nu, t);
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 8; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 8; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(15) << K(i, j);
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;
    
    return 0;
}