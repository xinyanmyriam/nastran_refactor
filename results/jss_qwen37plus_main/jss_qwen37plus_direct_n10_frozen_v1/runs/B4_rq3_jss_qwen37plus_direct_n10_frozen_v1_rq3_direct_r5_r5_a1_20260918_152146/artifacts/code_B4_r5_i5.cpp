#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

// Using double precision throughout
using Real = double;
using Matrix8d = Eigen::Matrix<Real, 8, 8>;
using Vector2d = Eigen::Vector2d;
using Matrix2d = Eigen::Matrix<Real, 2, 2>;
using Matrix6d = Eigen::Matrix<Real, 6, 6>;

// Helper: compute 2x2 rotation matrix for angle theta
Matrix2d rotationMatrix(Real theta) {
    return (Matrix2d() << std::cos(theta), -std::sin(theta),
                           std::sin(theta),  std::cos(theta)).finished();
}

// Compute full 6x6 stiffness for a CST triangle in GLOBAL coordinates
// But we will rotate it — so compute it in *local* frame defined by edge direction.
// Input: p0 = centroid, p1 = node[i], p2 = node[(i+1)%4]
// Local x-axis = (p2 - p1) normalized
Matrix6d computeTriangularStiffnessLocal(
    const Vector2d& p0, const Vector2d& p1, const Vector2d& p2,
    Real E, Real nu, Real t) {
    
    // Material matrix D for plane stress
    Real denom = 1.0 - nu * nu;
    Matrix2d D;
    D << E / denom, E * nu / denom,
         E * nu / denom, E / denom;

    // Area of triangle
    Real area = 0.5 * std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) -
                               (p2.x() - p0.x()) * (p1.y() - p0.y()));
    if (area <= 0.0) {
        return Matrix6d::Zero();
    }

    // Node coordinates relative to p0 (centroid) for shape function derivs
    Real x1 = p1.x() - p0.x(), y1 = p1.y() - p0.y();
    Real x2 = p2.x() - p0.x(), y2 = p2.y() - p0.y();

    Real detJ = x1 * y2 - x2 * y1;
    if (std::abs(detJ) < 1e-15) {
        return Matrix6d::Zero();
    }

    // b_i, c_i for shape functions N_i = a_i + b_i*x + c_i*y (in local coord centered at p0)
    // Standard CST: b1 = y2 - y1? Wait — re-derive for nodes (0,1,2) = (p0,p1,p2)
    // Let node0 = p0, node1 = p1, node2 = p2.
    // Then area = 0.5 |detJ|, and
    // b0 = y1 - y2, c0 = x2 - x1
    // b1 = y2 - y0 = y2, c1 = x0 - x2 = -x2
    // b2 = y0 - y1 = -y1, c2 = x1 - x0 = x1
    // But easier: use standard formula with (x0,y0), (x1,y1), (x2,y2)
    Real x0 = 0.0, y0 = 0.0; // p0 is origin in local frame
    Real b0 = y1 - y2;
    Real b1 = y2 - y0;
    Real b2 = y0 - y1;
    Real c0 = x2 - x1;
    Real c1 = x0 - x2;
    Real c2 = x1 - x0;

    Real inv2A = 1.0 / (2.0 * area);
    Real factor = t * area * inv2A * inv2A;

    Matrix6d K_local = Matrix6d::Zero();

    // Loop over node pairs i,j ∈ {0,1,2}
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Real b_i = (i == 0) ? b0 : (i == 1) ? b1 : b2;
            Real c_i = (i == 0) ? c0 : (i == 1) ? c1 : c2;
            Real b_j = (j == 0) ? b0 : (j == 1) ? b1 : b2;
            Real c_j = (j == 0) ? c0 : (j == 1) ? c1 : c2;

            // K(2*i, 2*j) += factor * (D00*b_i*b_j + D01*b_i*c_j + D10*c_i*b_j + D11*c_i*c_j)
            K_local(2*i, 2*j) += factor * (D(0,0)*b_i*b_j + D(0,1)*b_i*c_j + D(1,0)*c_i*b_j + D(1,1)*c_i*c_j);
            // K(2*i, 2*j+1) += factor * (D00*b_i*c_j + D01*b_i*b_j + D10*c_i*c_j + D11*c_i*b_j)
            K_local(2*i, 2*j+1) += factor * (D(0,0)*b_i*c_j + D(0,1)*b_i*b_j + D(1,0)*c_i*c_j + D(1,1)*c_i*b_j);
            // K(2*i+1, 2*j) += factor * (D00*c_i*b_j + D01*c_i*c_j + D10*b_i*b_j + D11*b_i*c_j)
            K_local(2*i+1, 2*j) += factor * (D(0,0)*c_i*b_j + D(0,1)*c_i*c_j + D(1,0)*b_i*b_j + D(1,1)*b_i*c_j);
            // K(2*i+1, 2*j+1) += factor * (D00*c_i*c_j + D01*c_i*b_j + D10*b_i*c_j + D11*b_i*b_j)
            K_local(2*i+1, 2*j+1) += factor * (D(0,0)*c_i*c_j + D(0,1)*c_i*b_j + D(1,0)*b_i*c_j + D(1,1)*b_i*b_j);
        }
    }

    return K_local;
}

// Assemble 8x8 stiffness matrix for CQUAD4 using 4 centroid-based triangles
Matrix8d computeCQUAD4Stiffness(
    const std::vector<Vector2d>& nodes, // 4 nodes: CCW order expected
    Real E, Real nu, Real t) {

    // Centroid
    Vector2d centroid = Vector2d::Zero();
    for (const auto& n : nodes) centroid += n;
    centroid /= 4.0;

    // Four triangles: centroid + nodes[i] + nodes[(i+1)%4], CCW
    // Ensure CCW winding: for convex quad, this is safe if input is CCW.
    std::vector<std::vector<int>> triangles = {
        {0, 1}, // tri0: centroid, nodes[0], nodes[1]
        {1, 2}, // tri1: centroid, nodes[1], nodes[2]
        {2, 3}, // tri2: centroid, nodes[2], nodes[3]
        {3, 0}  // tri3: centroid, nodes[3], nodes[0]
    };

    Matrix8d K = Matrix8d::Zero();

    for (int tri_idx = 0; tri_idx < 4; ++tri_idx) {
        int i = triangles[tri_idx][0];
        int j = triangles[tri_idx][1];
        const Vector2d& p1 = nodes[i];
        const Vector2d& p2 = nodes[j];

        // Local x-axis = direction from p1 to p2
        Vector2d edge = p2 - p1;
        Real len = edge.norm();
        if (len < 1e-15) continue;
        Vector2d ex = edge / len;
        Vector2d ey(-ex.y(), ex.x()); // CCW 90° rotation → local y

        // Angle of ex w.r.t global X
        Real theta = std::atan2(ex.y(), ex.x());
        Matrix2d T_node = rotationMatrix(theta); // maps (u,v) → (ux,uy)

        // Build 6x6 T = block-diag(T_node, T_node, T_node) for 3 nodes: centroid, p1, p2
        // Local DOF order: [u_c,v_c, u_p1,v_p1, u_p2,v_p2]
        // Global DOF order: centroid has no DOF! So we only map p1 and p2.
        // But centroid is not a physical node → its DOFs are *eliminated*.
        // So triangle has only 2 physical nodes → 4 DOFs, not 6.
        // Correction: In centroid-decomposition, centroid is *not assigned DOFs*. Only corner nodes have DOFs.
        // So each triangle contributes to exactly 2 corner nodes → 4 DOFs.
        // Therefore, local stiffness is 4x4, not 6x6.
        // We compute K_local for nodes p1 and p2 only (skip centroid).
        // So redefine triangle as 2-node edge + centroid → shape functions involve all 3, but only p1/p2 displace.

        // ✅ Correct approach: CST triangle with nodes [C, A, B] contributes to DOFs of A and B only.
        // So local DOF indices: 0→A, 1→B → local size = 4.
        // Compute 4x4 K_local for nodes A and B (indices 1 and 2 in triangle).
        // So we compute full 6x6 then extract rows/cols for local indices 1 and 2 (i.e., global nodes i and j).

        const Vector2d& p0 = centroid;
        Matrix6d K_local_full = computeTriangularStiffnessLocal(p0, p1, p2, E, nu, t);

        // Extract 4x4 submatrix for nodes 1 and 2 (local indices 1,2 → global i,j)
        // Local DOF order: [C, A, B] → DOFs: [cux,cuy, aux,auy, bux,buy]
        // We want rows/cols 2,3,4,5 → corresponding to A and B.
        Matrix6d K_local_4x4_full = Matrix6d::Zero();
        K_local_4x4_full.block<2,2>(0,0) = K_local_full.block<2,2>(2,2); // A-A
        K_local_4x4_full.block<2,2>(0,2) = K_local_full.block<2,2>(2,4); // A-B
        K_local_4x4_full.block<2,2>(2,0) = K_local_full.block<2,2>(4,2); // B-A
        K_local_4x4_full.block<2,2>(2,2) = K_local_full.block<2,2>(4,4); // B-B

        // Now build 4x4 T matrix: block-diag(T_node, T_node) for A and B
        Matrix2d T4 = Matrix2d::Zero();
        T4.block<2,2>(0,0) = T_node;
        T4.block<2,2>(2,2) = T_node;
        // Actually, T for 4 DOFs is:
        Eigen::Matrix<Real, 4, 4> T = Eigen::Matrix<Real, 4, 4>::Zero();
        T.block<2,2>(0,0) = T_node;
        T.block<2,2>(2,2) = T_node;

        // Rotate: K_global_4x4 = T^T * K_local_4x4 * T
        Eigen::Matrix<Real, 4, 4> K_global_4x4 = T.transpose() * K_local_4x4_full.block<4,4>(0,0) * T;

        // Map to global DOFs:
        // Node i → DOFs [2*i, 2*i+1]
        // Node j → DOFs [2*j, 2*j+1]
        int dofs_i[2] = {2*i, 2*i+1};
        int dofs_j[2] = {2*j, 2*j+1};

        for (int di = 0; di < 2; ++di) {
            for (int dj = 0; dj < 2; ++dj) {
                // A-A block: di,dj → [di,dj]
                K(dofs_i[di], dofs_i[dj]) += K_global_4x4(di, dj);
                // A-B block: di,dj → [di, 2+dj]
                K(dofs_i[di], dofs_j[dj]) += K_global_4x4(di, 2+dj);
                // B-A block: di,dj → [2+di, dj]
                K(dofs_j[di], dofs_i[dj]) += K_global_4x4(2+di, dj);
                // B-B block: di,dj → [2+di, 2+dj]
                K(dofs_j[di], dofs_j[dj]) += K_global_4x4(2+di, 2+dj);
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
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    return s;
}

int main() {
    // Test case inputs: CCW square-like quad
    std::vector<Vector2d> nodes = {
        Vector2d(0.0, 0.0),   // Node 0
        Vector2d(2.0, 0.0),   // Node 1
        Vector2d(2.0, 1.5),   // Node 2
        Vector2d(0.0, 1.5)    // Node 3
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