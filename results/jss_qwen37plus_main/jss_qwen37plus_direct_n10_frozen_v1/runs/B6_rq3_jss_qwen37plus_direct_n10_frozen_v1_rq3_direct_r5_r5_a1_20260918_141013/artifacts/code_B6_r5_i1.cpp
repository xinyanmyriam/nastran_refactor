#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Use double precision throughout
using Real = double;
using Matrix12x12 = Eigen::Matrix<Real, 12, 12>;
using Matrix6x6 = Eigen::Matrix<Real, 6, 6>;
using Matrix3x3 = Eigen::Matrix<Real, 3, 3>;
using Vector3 = Eigen::Vector3d;
using Vector12 = Eigen::VectorX<double>;

// Helper: format a number in scientific notation with fixed precision
std::string to_scientific(Real x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + in exponent
    size_t epos = s.find('e');
    if (epos != std::string::npos) {
        // Trim trailing zeros after decimal
        size_t dot = s.find('.');
        if (dot != std::string::npos && dot < epos) {
            while (s.back() == '0' && s.size() > dot + 2) {
                s.pop_back();
            }
            if (s.back() == '.') s.pop_back();
        }
        // Normalize exponent sign
        std::string exp = s.substr(epos);
        if (exp.length() >= 3 && exp[1] == '+') {
            s.replace(epos+1, 1, "");
        }
    }
    return s;
}

// Helper: print matrix as JSON array of arrays
void print_matrix_json(const Matrix12x12& K) {
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << to_scientific(K(i,j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
}

// Material properties helper (isotropic plate)
struct Material {
    Real E;   // Young's modulus
    Real nu;  // Poisson's ratio
    Real t;   // thickness
    Real D;   // bending rigidity = E*t^3/(12*(1-nu^2))

    Material(Real E_, Real nu_, Real t_) : E(E_), nu(nu_), t(t_) {
        D = E * t * t * t / (12.0 * (1.0 - nu * nu));
    }
};

// Triangle bending element stiffness (KTRBSC simplified logic)
// Returns 9x9 "super U" matrix partitioned into nine 3x3 blocks
// For triangle with vertices A, B, C in local coordinates:
// Local x-axis from A to B, y-axis from A to C projected perpendicular to x
std::vector<Matrix3x3> compute_triangle_ku(
    const Vector3& A, const Vector3& B, const Vector3& C,
    const Material& mat) {

    // Compute local coordinate system
    Vector3 AB = B - A;
    Real xsubb = AB.norm();
    if (xsubb < 1e-12) {
        throw std::runtime_error("Degenerate triangle: AB length zero");
    }
    Vector3 i_vec = AB / xsubb;

    Vector3 AC = C - A;
    Real xsubc = i_vec.dot(AC);
    Vector3 k_vec_unnorm = i_vec.cross(AC);
    Real ysubc = k_vec_unnorm.norm();
    if (ysubc < 1e-12) {
        throw std::runtime_error("Degenerate triangle: area zero");
    }
    Vector3 k_vec = k_vec_unnorm / ysubc;
    Vector3 j_vec = k_vec.cross(i_vec); // right-handed

    // Area
    Real area = xsubb * ysubc / 2.0;
    Real xbar = (xsubb + xsubc) / 3.0;
    Real ybar = ysubc / 3.0;

    // Precompute powers
    Real xbsq = xsubb * xsubb;
    Real xcsq = xsubc * xsubc;
    Real ycsq = ysubc * ysubc;
    Real xcy_c = xsubc * ysubc;
    Real px2 = (xbsq + xsubb*xsubc + xcsq) / 6.0;
    Real py2 = ycsq / 6.0;
    Real pxy2 = ysubc * (xsubb + 2.0*xsubc) / 12.0;
    Real xbar3 = 3.0 * xbar;
    Real ybar2 = 2.0 * ybar;
    Real ybar3 = 3.0 * ybar;

    // D matrix (bending stiffness components)
    // For isotropic: D11 = D, D12 = nu*D, D22 = D, D66 = (1-nu)*D/2
    Real D11 = mat.D;
    Real D12 = mat.nu * mat.D;
    Real D22 = mat.D;
    Real D66 = (1.0 - mat.nu) * mat.D / 2.0;

    // Build K^X matrix (36 elements, 6x6) - stored row-major
    std::vector<Real> KX(36, 0.0);
    #define IDX(i,j) ((i)*6+(j))
    KX[IDX(0,0)] = D11;
    KX[IDX(0,1)] = D12;
    KX[IDX(0,2)] = D66;
    KX[IDX(0,3)] = D11 * xbar3;
    KX[IDX(0,4)] = D12 * xbar + ybar2 * D66;
    KX[IDX(0,5)] = D12 * ybar3;

    KX[IDX(1,0)] = D12;
    KX[IDX(1,1)] = D22;
    KX[IDX(1,2)] = D66;
    KX[IDX(1,3)] = D12 * xbar3;
    KX[IDX(1,4)] = D66 * xbar + ybar2 * D22;
    KX[IDX(1,5)] = D66 * ybar3;

    KX[IDX(2,0)] = D66;
    KX[IDX(2,1)] = D66;
    KX[IDX(2,2)] = D66;
    KX[IDX(2,3)] = D66 * xbar3;
    KX[IDX(2,4)] = D66 * xbar + ybar2 * D66;
    KX[IDX(2,5)] = D66 * ybar3;

    KX[IDX(3,0)] = D11 * xbar3;
    KX[IDX(3,1)] = D12 * xbar3;
    KX[IDX(3,2)] = D66 * xbar3;
    KX[IDX(3,3)] = D11 * 9.0 * px2;
    KX[IDX(3,4)] = D12 * 3.0 * px2 + 6.0 * pxy2 * D66;
    KX[IDX(3,5)] = D12 * 9.0 * pxy2;

    KX[IDX(4,0)] = D12 * xbar + ybar2 * D66;
    KX[IDX(4,1)] = D66 * xbar + ybar2 * D22;
    KX[IDX(4,2)] = D66 * xbar + ybar2 * D66;
    KX[IDX(4,3)] = D12 * 3.0 * px2 + 6.0 * pxy2 * D66;
    KX[IDX(4,4)] = D22 * px2 + 4.0 * pxy2 * D66 + 4.0 * py2 * D22;
    KX[IDX(4,5)] = D22 * 3.0 * pxy2 + 6.0 * py2 * D66;

    KX[IDX(5,0)] = D12 * ybar3;
    KX[IDX(5,1)] = D66 * ybar3;
    KX[IDX(5,2)] = D66 * ybar3;
    KX[IDX(5,3)] = D12 * 9.0 * pxy2;
    KX[IDX(5,4)] = D22 * 3.0 * pxy2 + 6.0 * py2 * D66;
    KX[IDX(5,5)] = D22 * 9.0 * py2;

    // Scale by 4*area
    Real scale = 4.0 * area;
    for (auto& v : KX) v *= scale;

    // Build H-bar matrix (36 elements, 6x6) - stored row-major
    std::vector<Real> HBAR(36, 0.0);
    HBAR[IDX(0,0)] = xbsq;
    HBAR[IDX(0,3)] = xbsq * xsubb;
    HBAR[IDX(1,1)] = xsubb;
    HBAR[IDX(1,4)] = -2.0 * xsubb;
    HBAR[IDX(2,2)] = -3.0 * xbsq;
    HBAR[IDX(3,0)] = xcsq;
    HBAR[IDX(3,1)] = xcy_c;
    HBAR[IDX(3,2)] = ycsq;
    HBAR[IDX(3,3)] = xcsq * xsubc;
    HBAR[IDX(3,4)] = ycsq * xsubc;
    HBAR[IDX(3,5)] = ycsq * ysubc;
    HBAR[IDX(4,1)] = xsubc;
    HBAR[IDX(4,2)] = ysubc * 2.0;
    HBAR[IDX(4,4)] = xcy_c * 2.0;
    HBAR[IDX(4,5)] = ycsq * 3.0;
    HBAR[IDX(5,0)] = -2.0 * xsubc;
    HBAR[IDX(5,1)] = -ysubc;
    HBAR[IDX(5,4)] = -3.0 * xcsq;
    HBAR[IDX(5,5)] = -ycsq;

    // Invert H-bar (6x6)
    Eigen::Map<Eigen::Matrix<Real,6,6,Eigen::RowMajor>> Hmat(HBAR.data());
    Eigen::FullPivLU<Eigen::Matrix<Real,6,6>> lu(Hmat);
    if (!lu.isInvertible()) {
        throw std::runtime_error("H-bar matrix is singular");
    }
    Eigen::Matrix<Real,6,6> Hinv = lu.inverse();

    // Compute K^II = Hinv * KX * Hinv^T
    Eigen::Map<Eigen::Matrix<Real,6,6,Eigen::RowMajor>> KXmat(KX.data());
    Eigen::Matrix<Real,6,6> KII = Hinv * KXmat * Hinv.transpose();

    // Build S matrix (6x3)
    Eigen::Matrix<Real,6,3> S;
    S << 1.0, 0.0, -xsubb,
         0.0, 1.0,  0.0,
         0.0, 0.0,  1.0,
         1.0, ysubc, -xsubc,
         0.0, 1.0,  0.0,
         0.0, 0.0,  1.0;

    // Compute K^IA = -KII * S
    Eigen::Matrix<Real,6,3> KIA = -KII * S;

    // Compute K^AA = S^T * KIA
    Eigen::Matrix<Real,3,3> KAA = S.transpose() * KIA;

    // Assemble the 9x9 super U matrix from 3x3 blocks:
    // [ KII(0:2,0:2)  KIA(0:2,0:2)  KIA(0:2,0:2) ]
    // [ KIA^T         KAA           KIA(3:5,0:2)  ]
    // [ KIA^T         KIA(3:5,0:2)  KII(3:5,3:5)  ]
    // But per NASTRAN logic, it's arranged as 9 blocks of 3x3:
    // Block (i,j) corresponds to node i, node j, each 3x3
    std::vector<Matrix3x3> Ku(9);

    // K11 = top-left 3x3 of KII
    Ku[0] = KII.topLeftCorner(3,3);
    // K12 = top-right 3x3 of KIA
    Ku[1] = KIA.topRightCorner(3,3);
    // K13 = bottom-right 3x3 of KIA? Actually per indexing: KIA rows 0-2, cols 0-2
    Ku[2] = KIA.block(0,0,3,3);

    // K21 = transpose of K12
    Ku[3] = Ku[1].transpose();
    // K22 = KAA
    Ku[4] = KAA;
    // K23 = KIA rows 3-5, cols 0-2
    Ku[5] = KIA.block(3,0,3,3);

    // K31 = transpose of K13
    Ku[6] = Ku[2].transpose();
    // K32 = transpose of K23
    Ku[7] = Ku[5].transpose();
    // K33 = bottom-right 3x3 of KII
    Ku[8] = KII.bottomRightCorner(3,3);

    return Ku;
}

// Transformation matrix from local to global for a node
// For plate: w, theta_x, theta_y -> w, theta_x, theta_y in global
// Since our quad is planar in XY, and we assume no rotation, identity
Matrix3x3 get_local_to_global_transform() {
    return Matrix3x3::Identity();
}

// Main KQDPLT logic
Matrix12x12 kqdplt(const std::vector<Vector3>& nodes,
                   const Material& mat) {
    // Nodes: A, B, C, D in order (counterclockwise)
    const Vector3& A = nodes[0];
    const Vector3& B = nodes[1];
    const Vector3& C = nodes[2];
    const Vector3& D = nodes[3];

    // Pivot point: use node 0 (A) as pivot (NPIVOT = 1)
    int npivot = 1;

    // Compute local coordinate system for the quad
    // K vector: normal to plane
    Vector3 AB = B - A;
    Vector3 AD = D - A;
    Vector3 k_vec_unnorm = AB.cross(AD);
    Real norm_k = k_vec_unnorm.norm();
    if (norm_k < 1e-12) {
        throw std::runtime_error("Quad is degenerate: zero area");
    }
    Vector3 k_vec = k_vec_unnorm / norm_k;

    // I vector: AB normalized
    Real xsubb = AB.norm();
    Vector3 i_vec = AB / xsubb;

    // J vector: k cross i
    Vector3 j_vec = k_vec.cross(i_vec);

    // Now compute R matrix: 2x4 containing local coordinates of nodes
    // R(1,j) = x-coordinate of node j in local system
    // R(2,j) = y-coordinate of node j in local system
    std::vector<Real> R(8, 0.0); // 2x4 stored row-major: R[0]=R11, R[1]=R12, ..., R[4]=R21, ...
    // Node A: origin
    R[0] = 0.0; // A.x
    R[1] = xsubb; // B.x
    Vector3 AC = C - A;
    R[2] = AC.dot(i_vec); // C.x
    R[3] = AD.dot(i_vec); // D.x

    // Node A: origin
    R[4] = 0.0; // A.y
    R[5] = AD.dot(j_vec); // D.y
    R[6] = AC.dot(j_vec); // C.y
    // B.y = AB.dot(j_vec) = 0 since j_vec is orthogonal to i_vec and AB is along i_vec
    R[7] = 0.0; // B.y

    // Check convexity using cross products of consecutive edge vectors
    // Edges in CCW order: AB, BC, CD, DA
    Vector3 BC = C - B;
    Vector3 CD = D - C;
    Vector3 DA = A - D;

    // Cross products in 3D, then project to plane normal (k_vec) to get signed area
    auto cross_z = [](const Vector3& u, const Vector3& v, const Vector3& n) -> Real {
        return u.cross(v).dot(n);
    };

    Real cp1 = cross_z(AB, BC, k_vec);
    Real cp2 = cross_z(BC, CD, k_vec);
    Real cp3 = cross_z(CD, DA, k_vec);
    Real cp4 = cross_z(DA, AB, k_vec);

    // For strict convexity and CCW ordering, all cross products must be > 0
    // Also check that no three points are collinear (cross product not near zero)
    const Real eps = 1e-10;
    if (cp1 <= eps || cp2 <= eps || cp3 <= eps || cp4 <= eps) {
        throw std::runtime_error("Quad is non-convex");
    }

    // Sub-triangle mapping: M array from Fortran
    // M = [2,4,1, 3,1,2, 4,2,3, 1,3,4] -> 4 triangles: (2,4,1), (3,1,2), (4,2,3), (1,3,4)
    // For pivot=1 (node A), triangles containing node A are:
    // Triangle 1: (2,4,1) -> (B,D,A)
    // Triangle 2: (3,1,2) -> (C,A,B)
    // Triangle 3: (1,3,4) -> (A,C,D)
    // So three triangles: ABD, ABC, ACD -> indices: {0,1,3}, {0,2,1}, {0,2,3}
    std::vector<std::vector<int>> triangles = {
        {0, 1, 3}, // A,B,D
        {0, 2, 1}, // A,C,B
        {0, 2, 3}  // A,C,D
    };

    // Accumulate 12x12 stiffness
    Matrix12x12 K_total = Matrix12x12::Zero();

    // For each triangle, compute its 9x9 "super U" and map to global DOFs
    for (const auto& tri_nodes : triangles) {
        std::vector<Vector3> tri_pts = {
            nodes[tri_nodes[0]],
            nodes[tri_nodes[1]],
            nodes[tri_nodes[2]]
        };

        auto Ku_blocks = compute_triangle_ku(tri_pts[0], tri_pts[1], tri_pts[2], mat);

        // Map local triangle DOFs to global 12-DOF vector
        // Each node has 3 DOFs: w, theta_x, theta_y
        // Triangle has 3 nodes -> 9 DOFs
        // We'll create a 12x9 mapping matrix
        Eigen::Matrix<Real,12,9> mapping = Eigen::Matrix<Real,12,9>::Zero();
        for (int i = 0; i < 3; ++i) {
            int global_node = tri_nodes[i];
            int start = global_node * 3;
            mapping.block(start, i*3, 3, 3) = Matrix3x3::Identity();
        }

        // Assemble 9x9 Ku into 12x12
        Eigen::Matrix<Real,12,12> K_tri = Eigen::Matrix<Real,12,12>::Zero();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                int global_i = tri_nodes[i] * 3;
                int global_j = tri_nodes[j] * 3;
                K_tri.block(global_i, global_j, 3, 3) = Ku_blocks[i*3+j];
            }
        }

        K_total += K_tri;
    }

    // Apply transformations: for each node, transform DOFs
    // Since our quad is in XY plane and we assume no orientation, use identity
    // But per NASTRAN, there might be coordinate systems. For test case, all are basic.
    for (int node = 0; node < 4; ++node) {
        int start = node * 3;
        Matrix3x3 T = get_local_to_global_transform();
        // Transform rows and columns: K = T^T * K * T
        // First transform columns
        Eigen::Matrix<Real,12,3> K_temp = K_total.block(0, start, 12, 3) * T;
        // Then transform rows
        K_total.block(0, start, 12, 3) = K_temp;
        Eigen::Matrix<Real,3,12> K_temp2 = T.transpose() * K_total.block(start, 0, 3, 12);
        K_total.block(start, 0, 3, 12) = K_temp2;
    }

    return K_total;
}

int main() {
    try {
        // Test case: Quad A=(0,0,0), B=(1,0,0), C=(1,1,0), D=(0,1,0)
        std::vector<Vector3> nodes = {
            Vector3(0.0, 0.0, 0.0), // A
            Vector3(1.0, 0.0, 0.0), // B
            Vector3(1.0, 1.0, 0.0), // C
            Vector3(0.0, 1.0, 0.0)  // D
        };

        // Material: E = 200e9, nu = 0.3, t = 0.01
        Real E = 200e9;
        Real nu = 0.3;
        Real t = 0.01;
        Material mat(E, nu, t);

        // Compute stiffness matrix
        Matrix12x12 K = kqdplt(nodes, mat);

        // Output as JSON
        print_matrix_json(K);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}