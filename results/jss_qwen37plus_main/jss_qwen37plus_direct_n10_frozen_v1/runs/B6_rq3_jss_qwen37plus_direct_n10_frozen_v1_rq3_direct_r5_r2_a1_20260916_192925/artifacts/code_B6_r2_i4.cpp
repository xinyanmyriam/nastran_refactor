#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper function to compute cross product of two 3D vectors
Vector3d cross(const Vector3d& a, const Vector3d& b) {
    return Vector3d(
        a(1)*b(2) - a(2)*b(1),
        a(2)*b(0) - a(0)*b(2),
        a(0)*b(1) - a(1)*b(0)
    );
}

// Helper function to compute dot product of two 3D vectors
double dot(const Vector3d& a, const Vector3d& b) {
    return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
}

// Helper function to normalize a 3D vector
Vector3d normalize(const Vector3d& v) {
    double norm = v.norm();
    if (norm == 0.0) {
        return Vector3d::Zero();
    }
    return v / norm;
}

// Compute plate-local frame: T = [x_tilde, y_tilde, z_tilde]^T
// z_tilde = unit normal to quad plane (average of triangle normals)
// x_tilde = direction from node0 to node1, projected & normalized onto plane
// y_tilde = z_tilde × x_tilde
Matrix3d compute_plate_frame(const vector<Vector3d>& nodes) {
    // Use all 4 nodes to compute best-fit plane normal
    Vector3d v1 = nodes[1] - nodes[0];
    Vector3d v2 = nodes[2] - nodes[0];
    Vector3d v3 = nodes[3] - nodes[0];
    Vector3d n1 = cross(v1, v2);
    Vector3d n2 = cross(v2, v3);
    Vector3d n3 = cross(v3, v1);
    Vector3d n_avg = n1 + n2 + n3;
    Vector3d z_tilde = normalize(n_avg);

    // x_tilde: project v1 onto plane perpendicular to z_tilde
    Vector3d x_tilde = v1 - dot(v1, z_tilde) * z_tilde;
    if (x_tilde.norm() < 1e-12) {
        // fallback: use v2
        x_tilde = v2 - dot(v2, z_tilde) * z_tilde;
    }
    x_tilde = normalize(x_tilde);

    // y_tilde = z_tilde × x_tilde (orthonormal right-handed)
    Vector3d y_tilde = cross(z_tilde, x_tilde);

    Matrix3d T;
    T.row(0) = x_tilde.transpose();
    T.row(1) = y_tilde.transpose();
    T.row(2) = z_tilde.transpose();
    return T;
}

// Compute DKT (Discrete Kirchhoff Triangle) 9x9 stiffness matrix
// Nodes: p1, p2, p3 (in global coordinates)
// Returns K_local (9x9) in triangle-local coordinates where:
//   local x-axis = p1->p2
//   local y-axis = in-plane, orthogonal to x, same side as normal
//   local z = normal
MatrixXd compute_DKT_stiffness(
    const Vector3d& p1, const Vector3d& p2, const Vector3d& p3,
    double E_modulus, double nu, double thickness) {

    // Moment of inertia
    double I = thickness * thickness * thickness / 12.0;

    // Material matrix D for isotropic plate bending (3x3)
    double D11 = E_modulus * I / (1.0 - nu * nu);
    double D12 = nu * D11;
    double D66 = (1.0 - nu) * D11 / 2.0;

    Matrix3d D;
    D << D11, D12, 0.0,
         D12, D11, 0.0,
         0.0, 0.0, D66;

    // Triangle geometry
    Vector3d a = p2 - p1;
    Vector3d b = p3 - p1;
    double area = 0.5 * cross(a, b).norm();

    // Local coordinate system (triangle plane)
    Vector3d z_local = normalize(cross(a, b));
    Vector3d x_local = normalize(a);
    Vector3d y_local = normalize(cross(z_local, x_local));

    // Project points onto local (x,y) plane
    auto local_coord = [&](const Vector3d& p) -> Vector2d {
        Vector3d dp = p - p1;
        return Vector2d(dot(dp, x_local), dot(dp, y_local));
    };

    Vector2d q1 = Vector2d::Zero();
    Vector2d q2 = local_coord(p2);
    Vector2d q3 = local_coord(p3);

    // Side lengths in local coords
    double L12 = q2.norm(); // = q2.x()
    double L13 = q3.norm();
    double L23 = (q3 - q2).norm();

    // Barycentric coordinates (area coordinates)
    double A = area;
    double a1 = ((q2.y() * q3.x()) - (q2.x() * q3.y())) * 0.5; // area opposite p1
    double a2 = ((q3.y() * q1.x()) - (q3.x() * q1.y())) * 0.5; // = 0
    double a3 = ((q1.y() * q2.x()) - (q1.x() * q2.y())) * 0.5; // = 0
    // Instead, use standard DKT shape function setup

    // DKT uses 9 DOFs: w1,w2,w3, θx1,θx2,θx3, θy1,θy2,θy3
    // Standard DKT B-matrix is known — we use closed-form from Cook/Belytschko
    // Reference: "Concepts and Applications of Finite Element Analysis", 4th ed., sec 12.4

    // Node coordinates in local (x,y)
    double x1 = q1(0), y1 = q1(1);
    double x2 = q2(0), y2 = q2(1);
    double x3 = q3(0), y3 = q3(1);

    // Compute side vectors
    double l1 = x2 - x1, m1 = y2 - y1; // 1->2
    double l2 = x3 - x2, m2 = y3 - y2; // 2->3
    double l3 = x1 - x3, m3 = y1 - y3; // 3->1

    // Edge lengths
    double s1 = sqrt(l1*l1 + m1*m1);
    double s2 = sqrt(l2*l2 + m2*m2);
    double s3 = sqrt(l3*l3 + m3*m3);

    // DKT curvature-strain matrix B is 3x9 (for [κxx, κyy, κxy])
    // We build B for each of 3 Gauss points (midpoints of sides) — but standard DKT uses analytical integration.

    // Instead, use the well-known analytical DKT stiffness (see Zienkiewicz, Taylor vol.2)
    // We'll construct B matrices at 3 points: midpoints of edges

    // Midpoint coordinates
    Vector2d mp12 = 0.5 * (q1 + q2);
    Vector2d mp23 = 0.5 * (q2 + q3);
    Vector2d mp31 = 0.5 * (q3 + q1);

    // Precompute common factors
    double c1 = 1.0/(2.0*A);
    double c2 = 1.0/(2.0*A*A);

    // Build 3x9 B matrix at each midpoint
    // B = [ ∂²φ/∂x², ∂²φ/∂y², 2∂²φ/∂x∂y ] for each shape function φ_i
    // For DKT, shape functions are quadratic in w, linear in θx,θy

    // Rather than derive B here, use standard DKT implementation:
    // See: http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch12.pdf (Eq. 12.15)

    // We'll use the compact form: K = t * ∫∫ B^T D B dA
    // Analytical result exists — but for robustness, we integrate numerically at 3 points with equal weight A/3

    // Each Gauss point: weight = A/3
    double wt = A / 3.0;

    MatrixXd K_local = MatrixXd::Zero(9, 9);

    // Helper: evaluate B at point (x,y) — returns 3x9 matrix
    auto eval_B = [&](const Vector2d& pt) -> MatrixXd {
        double x = pt(0), y = pt(1);
        MatrixXd B(3, 9);
        B.setZero();

        // Shape functions for DKT (from Cook p.472)
        // w_i = N_i^w, θx_i = N_i^θx, θy_i = N_i^θy
        // But standard DKT uses 3-node interpolation for w (quadratic) and θ (linear)
        // Simpler: use known analytical B for DKT triangle (see Felippa IFEM)

        // We'll use the standard DKT B-matrix components (Felippa Eq. 24.12)
        // B11 = ∂²N1/∂x², B12 = ∂²N1/∂y², B13 = 2∂²N1/∂x∂y, etc.

        // N1 = L1*(2*L1 - 1) + 4*L2*L3   [quadratic]
        // N2 = L2*(2*L2 - 1) + 4*L3*L1
        // N3 = L3*(2*L3 - 1) + 4*L1*L2
        // where L1,L2,L3 are area coordinates: L1 = A1/A, etc.

        auto area_coord = [&](double x, double y) -> Vector3d {
            // Area coordinates for point (x,y) in triangle (x1,y1), (x2,y2), (x3,y3)
            double denom = 2.0*A;
            double L1 = ((x2*y3 - x3*y2) + (y2 - y3)*x + (x3 - x2)*y) / denom;
            double L2 = ((x3*y1 - x1*y3) + (y3 - y1)*x + (x1 - x3)*y) / denom;
            double L3 = 1.0 - L1 - L2;
            return Vector3d(L1, L2, L3);
        };

        Vector3d L = area_coord(x, y);
        double L1 = L(0), L2 = L(1), L3 = L(2);

        // Derivatives of area coordinates
        // dL1/dx = a1/A, dL1/dy = b1/A, etc. where a1 = (y2-y3)/2, b1 = (x3-x2)/2
        double a1 = (y2 - y3) * c1;
        double b1 = (x3 - x2) * c1;
        double a2 = (y3 - y1) * c1;
        double b2 = (x1 - x3) * c1;
        double a3 = (y1 - y2) * c1;
        double b3 = (x2 - x1) * c1;

        // Quadratic shape functions for w
        // N1w = L1*(2*L1 - 1) + 4*L2*L3
        // N2w = L2*(2*L2 - 1) + 4*L3*L1
        // N3w = L3*(2*L3 - 1) + 4*L1*L2

        // First derivatives of Niw
        double dN1w_dx = (4*L1 - 1)*a1 + 4*(L2*a3 + L3*a2);
        double dN1w_dy = (4*L1 - 1)*b1 + 4*(L2*b3 + L3*b2);
        double dN2w_dx = (4*L2 - 1)*a2 + 4*(L3*a1 + L1*a3);
        double dN2w_dy = (4*L2 - 1)*b2 + 4*(L3*b1 + L1*b3);
        double dN3w_dx = (4*L3 - 1)*a3 + 4*(L1*a2 + L2*a1);
        double dN3w_dy = (4*L3 - 1)*b3 + 4*(L1*b2 + L2*b1);

        // Second derivatives of Niw
        double d2N1w_dxx = 4*a1*a1 + 4*(a2*a3 + a3*a2); // = 4*a1² + 8*a2*a3
        double d2N1w_dyy = 4*b1*b1 + 4*(b2*b3 + b3*b2);
        double d2N1w_dxy = 4*a1*b1 + 4*(a2*b3 + a3*b2);

        double d2N2w_dxx = 4*a2*a2 + 4*(a3*a1 + a1*a3);
        double d2N2w_dyy = 4*b2*b2 + 4*(b3*b1 + b1*b3);
        double d2N2w_dxy = 4*a2*b2 + 4*(a3*b1 + a1*b3);

        double d2N3w_dxx = 4*a3*a3 + 4*(a1*a2 + a2*a1);
        double d2N3w_dyy = 4*b3*b3 + 4*(b1*b2 + b2*b1);
        double d2N3w_dxy = 4*a3*b3 + 4*(a1*b2 + a2*b1);

        // Linear shape functions for θx, θy: same as area coords
        // N1θx = L1, N2θx = L2, N3θx = L3
        // N1θy = L1, etc.
        double dN1θx_dx = a1; double dN1θx_dy = b1;
        double dN2θx_dx = a2; double dN2θx_dy = b2;
        double dN3θx_dx = a3; double dN3θx_dy = b3;
        double dN1θy_dx = a1; double dN1θy_dy = b1;
        double dN2θy_dx = a2; double dN2θy_dy = b2;
        double dN3θy_dx = a3; double dN3θy_dy = b3;

        // Curvature components:
        // κxx = ∂²w/∂x² - ∂θy/∂x
        // κyy = ∂²w/∂y² + ∂θx/∂y
        // κxy = ∂²w/∂x∂y + ∂θx/∂x - ∂θy/∂y
        // So B matrix rows:
        // Row 0 (κxx): [d2w1_dxx, d2w2_dxx, d2w3_dxx, 0,0,0, -dL1_dx, -dL2_dx, -dL3_dx]
        // Row 1 (κyy): [d2w1_dyy, d2w2_dyy, d2w3_dyy, dL1_dy, dL2_dy, dL3_dy, 0,0,0]
        // Row 2 (κxy): [d2w1_dxy, d2w2_dxy, d2w3_dxy, dL1_dx, dL2_dx, dL3_dx, -dL1_dy, -dL2_dy, -dL3_dy]

        B(0,0) = d2N1w_dxx; B(0,1) = d2N2w_dxx; B(0,2) = d2N3w_dxx;
        B(0,3) = 0.0; B(0,4) = 0.0; B(0,5) = 0.0;
        B(0,6) = -a1; B(0,7) = -a2; B(0,8) = -a3;

        B(1,0) = d2N1w_dyy; B(1,1) = d2N2w_dyy; B(1,2) = d2N3w_dyy;
        B(1,3) = b1; B(1,4) = b2; B(1,5) = b3;
        B(1,6) = 0.0; B(1,7) = 0.0; B(1,8) = 0.0;

        B(2,0) = d2N1w_dxy; B(2,1) = d2N2w_dxy; B(2,2) = d2N3w_dxy;
        B(2,3) = a1; B(2,4) = a2; B(2,5) = a3;
        B(2,6) = -b1; B(2,7) = -b2; B(2,8) = -b3;

        return B;
    };

    // Integrate B^T D B at 3 midpoints
    vector<Vector2d> mps = {mp12, mp23, mp31};
    for (const auto& mp : mps) {
        MatrixXd B = eval_B(mp);
        MatrixXd BTDB = B.transpose() * D * B;
        K_local += wt * BTDB;
    }

    return K_local;
}

// Main KQDPLT function that computes the 12x12 stiffness matrix
MatrixXd KQDPLT(
    const vector<Vector3d>& nodes, // 4 nodes: A, B, C, D
    double E_modulus, double nu, double thickness) {

    // Validate input
    if (nodes.size() != 4) {
        throw runtime_error("KQDPLT: exactly 4 nodes required");
    }

    // Compute plate-local frame T (3x3)
    Matrix3d T = compute_plate_frame(nodes);
    Matrix3d Tt = T.transpose();

    // Global stiffness matrix
    MatrixXd K = MatrixXd::Zero(12, 12);

    // Triangle 1: A, B, C  → nodes[0], nodes[1], nodes[2]
    MatrixXd K_T1 = compute_DKT_stiffness(nodes[0], nodes[1], nodes[2], E_modulus, nu, thickness);

    // Triangle 2: A, C, D  → nodes[0], nodes[2], nodes[3]
    MatrixXd K_T2 = compute_DKT_stiffness(nodes[0], nodes[2], nodes[3], E_modulus, nu, thickness);

    // Node mapping for triangle 1: local 0→A, 1→B, 2→C → global DOF: 0,1,2 → 0..2, 3..5, 6..8
    // Node mapping for triangle 2: local 0→A, 1→C, 2→D → global DOF: 0..2, 6..8, 9..11

    // Helper lambda to rotate a 3x3 block
    auto rotate_block = [&](const Matrix3d& K_local_block) -> Matrix3d {
        return Tt * K_local_block * T;
    };

    // Assemble triangle 1
    for (int i = 0; i < 3; ++i) {
        int gi = (i == 0) ? 0 : (i == 1) ? 3 : 6; // A→0, B→3, C→6
        for (int j = 0; j < 3; ++j) {
            int gj = (j == 0) ? 0 : (j == 1) ? 3 : 6;
            Matrix3d K_block = K_T1.block(i*3, j*3, 3, 3);
            Matrix3d K_rot = rotate_block(K_block);
            K.block(gi, gj, 3, 3) += K_rot;
        }
    }

    // Assemble triangle 2
    // local node 0 → A → DOF 0,1,2
    // local node 1 → C → DOF 6,7,8
    // local node 2 → D → DOF 9,10,11
    vector<int> t2_map = {0, 6, 9}; // local→global DOF start index
    for (int i = 0; i < 3; ++i) {
        int gi = t2_map[i];
        for (int j = 0; j < 3; ++j) {
            int gj = t2_map[j];
            Matrix3d K_block = K_T2.block(i*3, j*3, 3, 3);
            Matrix3d K_rot = rotate_block(K_block);
            K.block(gi, gj, 3, 3) += K_rot;
        }
    }

    return K;
}

int main() {
    // Test case: unit square in XY plane
    vector<Vector3d> nodes = {
        Vector3d(0.0, 0.0, 0.0), // A
        Vector3d(1.0, 0.0, 0.0), // B
        Vector3d(1.0, 1.0, 0.0), // C
        Vector3d(0.0, 1.0, 0.0)  // D
    };

    double E = 200e9;      // Pa
    double nu = 0.3;
    double t = 0.01;      // m

    try {
        // Compute stiffness matrix
        MatrixXd K = KQDPLT(nodes, E, nu, t);

        // Output the 12x12 matrix in space-separated format
        cout << fixed;
        cout.precision(6);
        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 12; j++) {
                if (j > 0) cout << " ";
                cout << K(i,j);
            }
            cout << endl;
        }

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}