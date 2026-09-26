#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <stdexcept>
#include <sstream>

// For M_PI compatibility
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace Eigen;

// Constants
const double DEGRA = M_PI / 180.0;

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
        throw runtime_error("Cannot normalize zero vector");
    }
    return v / norm;
}

// Matrix multiplication helper: C = A * B^T (if transpose_b is true)
MatrixXd matmul(const MatrixXd& A, const MatrixXd& B, bool transpose_b = false) {
    if (transpose_b) {
        return A * B.transpose();
    } else {
        return A * B;
    }
}

// Matrix multiplication with transpose option for each matrix
MatrixXd gmmatd(const MatrixXd& A, int rows_a, int cols_a, int trans_a,
                const MatrixXd& B, int rows_b, int cols_b, int trans_b,
                bool transpose_result = false) {
    MatrixXd A_use = A;
    MatrixXd B_use = B;
    
    if (trans_a) {
        A_use = A.transpose();
    }
    if (trans_b) {
        B_use = B.transpose();
    }
    
    MatrixXd result = A_use * B_use;
    
    if (transpose_result) {
        result = result.transpose();
    }
    
    return result;
}

// Transformation matrix for coordinate system rotation
MatrixXd transd(int csid, const Vector3d& u1_vec, const Vector3d& u2_vec, const Vector3d& k_vec) {
    // For this implementation, we assume csid == 0 means no transformation
    // Otherwise, we use the provided basis vectors to form the transformation matrix
    if (csid == 0) {
        return MatrixXd::Identity(3, 3);
    }
    
    // Create transformation matrix from local to global coordinates
    // Columns are the local basis vectors expressed in global coordinates
    MatrixXd T(3, 3);
    T.col(0) = u1_vec;  // i-vector
    T.col(1) = u2_vec;  // j-vector  
    T.col(2) = k_vec;   // k-vector
    
    return T;
}

// Material properties routine (simplified for isotropic material)
void mat(double youngs_modulus, double nu, double& g11, double& g12, double& g13,
         double& g22, double& g23, double& g33) {
    // For isotropic material, the stiffness matrix components
    // G = [D] matrix for plate bending (reduced stiffness)
    double D = youngs_modulus / (1.0 - nu*nu);
    g11 = D;
    g12 = D * nu;
    g13 = 0.0;
    g22 = D;
    g23 = 0.0;
    g33 = D * (1.0 - nu) / 2.0;
}

// Basic bending triangle routine (KTRBSC)
void ktrbsc(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3,
            double youngs_modulus, double nu, double I, MatrixXd& K_super_u) {
    // Compute vectors
    Vector3d r_ab = p2 - p1;
    Vector3d r_ac = p3 - p1;
    
    // XSUBB = |r_ab|
    double xsubb = r_ab.norm();
    if (xsubb < 1e-12) {
        throw runtime_error("Triangle has zero area (degenerate)");
    }
    
    // I-vector = r_ab normalized
    Vector3d i_vec = r_ab / xsubb;
    
    // XSUBC = i_vec . r_ac
    double xsubc = dot(i_vec, r_ac);
    
    // K-vector = i_vec × r_ac (non-normalized)
    Vector3d k_vec = cross(i_vec, r_ac);
    double ysubc = k_vec.norm();
    if (ysubc < 1e-12) {
        throw runtime_error("Triangle has zero area (collinear points)");
    }
    
    // Normalize k-vector
    k_vec = k_vec / ysubc;
    
    // J-vector = k_vec × i_vec
    Vector3d j_vec = cross(k_vec, i_vec);
    j_vec = normalize(j_vec);
    
    // Material properties
    double g11, g12, g13, g22, g23, g33;
    mat(youngs_modulus, nu, g11, g12, g13, g22, g23, g33);
    
    // D matrix = I * G matrix
    Matrix3d D;
    D << g11*I, g12*I, g13*I,
         g12*I, g22*I, g23*I,
         g13*I, g23*I, g33*I;
    
    // Triangle area
    double area = xsubb * ysubc / 2.0;
    
    // Centroid coordinates in local system
    double xbar = (xsubb + xsubc) / 3.0;
    double ybar = ysubc / 3.0;
    
    // Precompute terms
    double xcsq = xsubc * xsubc;
    double ycsq = ysubc * ysubc;
    double xbsq = xsubb * xsubb;
    double xcy_c = xsubc * ysubc;
    double px2 = (xbsq + xsubb*xsubc + xcsq) / 6.0;
    double py2 = ycsq / 6.0;
    double pxy2 = ysubc * (xsubb + 2.0*xsubc) / 12.0;
    double xbar3 = 3.0 * xbar;
    double ybar2 = 2.0 * ybar;
    double ybar3 = 3.0 * ybar;
    
    // Build K^X matrix (36 elements, 6x6)
    MatrixXd KX(6, 6);
    KX << D(0,0), D(0,2), D(0,1), D(0,0)*xbar3, D(0,1)*xbar + ybar2*D(0,2), D(0,1)*ybar3,
           D(0,2), D(2,2), D(1,2), D(0,2)*xbar3, D(1,2)*xbar + ybar2*D(2,2), D(1,2)*ybar3,
           D(0,1), D(1,2), D(1,1), D(0,1)*xbar3, D(1,1)*xbar + ybar2*D(1,2), D(1,1)*ybar3,
           D(0,0)*xbar3, D(0,2)*xbar3, D(0,1)*xbar3, D(0,0)*9.0*px2, D(0,1)*3.0*px2 + 6.0*pxy2*D(0,2), D(0,1)*9.0*pxy2,
           D(0,1)*xbar + ybar2*D(0,2), D(1,2)*xbar + ybar2*D(2,2), D(1,1)*xbar + ybar2*D(1,2), D(0,1)*3.0*px2 + 6.0*pxy2*D(0,2), D(1,1)*px2 + 4.0*pxy2*D(1,2) + 4.0*py2*D(2,2), D(1,1)*3.0*pxy2 + 6.0*py2*D(1,2),
           D(0,1)*ybar3, D(1,2)*ybar3, D(1,1)*ybar3, D(0,1)*9.0*pxy2, D(1,1)*3.0*pxy2 + 6.0*py2*D(1,2), D(1,1)*9.0*py2;
    
    // Scale by 4*area
    KX *= 4.0 * area;
    
    // Build HBAR matrix (36 elements, 6x6)
    MatrixXd HBAR(6, 6);
    HBAR.setZero();
    HBAR(0,0) = xbsq;
    HBAR(1,0) = xbsq*xsubb;
    HBAR(3,0) = xsubb;
    HBAR(4,0) = -2.0*xsubb;
    HBAR(5,0) = -3.0*xbsq;
    HBAR(0,1) = xcsq;
    HBAR(1,1) = xcy_c;
    HBAR(2,1) = ycsq;
    HBAR(3,1) = xcsq*xsubc;
    HBAR(4,1) = ycsq*xsubc;
    HBAR(5,1) = ycsq*ysubc;
    HBAR(1,2) = xsubc;
    HBAR(2,2) = ysubc*2.0;
    HBAR(4,2) = xcy_c*2.0;
    HBAR(5,2) = ycsq*3.0;
    HBAR(0,3) = -2.0*xsubc;
    HBAR(1,3) = -ysubc;
    HBAR(5,3) = -3.0*xcsq;
    HBAR(2,3) = -ycsq;
    
    // Invert HBAR matrix
    MatrixXd HBAR_inv;
    try {
        HBAR_inv = HBAR.inverse();
    } catch (...) {
        throw runtime_error("HBAR matrix is singular");
    }
    
    // Compute KII = KX * HBAR_inv
    MatrixXd KII = KX * HBAR_inv;
    
    // Build S matrix (6x3)
    MatrixXd S(6, 3);
    S << 1.0, 0.0, -xsubb,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, ysubc, -xsubc,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    
    // Compute KIA = KII * S
    MatrixXd KIA = KII * S;
    
    // Compute KAA = S.transpose() * KIA
    MatrixXd KAA = S.transpose() * KIA;
    
    // Now arrange the 9 3x3 matrices of K super U
    // The structure is:
    // [KAA, KIA^T, KIA]
    // [KIA, KII, KII*S]
    // [KIA, S^T*KII, KII]
    
    // For simplicity, we'll construct K_super_u as a 9x9 matrix
    // But the original code arranges it as 9 3x3 blocks in a specific order
    // We'll create the 9x9 matrix directly
    
    K_super_u = MatrixXd::Zero(9, 9);
    
    // Fill KAA (top-left 3x3)
    K_super_u.block(0, 0, 3, 3) = KAA;
    
    // Fill KIA^T (top-middle 3x3)
    K_super_u.block(0, 3, 3, 3) = KIA.transpose();
    
    // Fill KIA (top-right 3x3)
    K_super_u.block(0, 6, 3, 3) = KIA;
    
    // Fill KIA (middle-left 3x3)
    K_super_u.block(3, 0, 3, 3) = KIA;
    
    // Fill KII (middle-middle 3x3)
    K_super_u.block(3, 3, 3, 3) = KII;
    
    // Fill KII*S (middle-right 3x3)
    K_super_u.block(3, 6, 3, 3) = KII * S;
    
    // Fill KIA (bottom-left 3x3)
    K_super_u.block(6, 0, 3, 3) = KIA;
    
    // Fill S^T*KII (bottom-middle 3x3)
    K_super_u.block(6, 3, 3, 3) = S.transpose() * KII;
    
    // Fill KII (bottom-right 3x3)
    K_super_u.block(6, 6, 3, 3) = KII;
}

// Main quadrilateral plate bending element routine (KQDPLT)
MatrixXd kqdplt(const vector<Vector3d>& nodes, double youngs_modulus, double nu, double t) {
    // Calculate moment of inertia
    double I = t*t*t / 12.0;
    
    // Node coordinates (A, B, C, D)
    Vector3d A = nodes[0];
    Vector3d B = nodes[1];
    Vector3d C = nodes[2];
    Vector3d D = nodes[3];
    
    // Determine pivot point (we'll use node A as pivot, index 1 in Fortran 1-based)
    int npivot = 1; // A is pivot (1-based indexing)
    
    // Check if pivot is valid
    if (npivot < 1 || npivot > 4) {
        throw runtime_error("Invalid pivot point");
    }
    
    // Compute normal vector using cross product of two diagonals or adjacent edges
    // Use (B-A) and (D-A) to get initial normal
    Vector3d k_vec = cross(B-A, D-A);
    double norm_k = k_vec.norm();
    if (norm_k < 1e-12) {
        throw runtime_error("Quadrilateral is degenerate (coplanar points not forming area)");
    }
    k_vec = k_vec / norm_k;
    
    // Construct orthonormal basis: i_vec along (B-A) projected to plane, j_vec = k_vec × i_vec
    Vector3d a1 = B - A;
    double proj_k = dot(a1, k_vec);
    Vector3d i_vec = a1 - proj_k * k_vec;
    double norm_i = i_vec.norm();
    if (norm_i < 1e-12) {
        // If B-A is parallel to k_vec, use another edge
        Vector3d a2 = C - A;
        proj_k = dot(a2, k_vec);
        i_vec = a2 - proj_k * k_vec;
        norm_i = i_vec.norm();
        if (norm_i < 1e-12) {
            throw runtime_error("Unable to construct local basis: all edges parallel to normal");
        }
    }
    i_vec = i_vec / norm_i;
    
    Vector3d j_vec = cross(k_vec, i_vec);
    j_vec = normalize(j_vec);
    
    // Compute coordinates in local system (x along i_vec, y along j_vec)
    // R(0,j) = x-coordinate, R(1,j) = y-coordinate
    MatrixXd R(2, 4);
    R(0,0) = 0.0; R(1,0) = 0.0; // A
    R(0,1) = dot(B-A, i_vec); R(1,1) = dot(B-A, j_vec); // B
    R(0,2) = dot(C-A, i_vec); R(1,2) = dot(C-A, j_vec); // C
    R(0,3) = dot(D-A, i_vec); R(1,3) = dot(D-A, j_vec); // D
    
    // Robust convexity check using 2D cross products
    // Points in order A,B,C,D (0,1,2,3). Compute cross product of consecutive edge vectors.
    // Edge vectors: AB, BC, CD, DA
    // Cross product z_i = (P[i+1]-P[i]) × (P[i+2]-P[i+1]) for i=0..3 (mod 4)
    vector<double> z(4);
    auto cross2d = [](const Vector2d& a, const Vector2d& b) -> double {
        return a(0)*b(1) - a(1)*b(0);
    };
    
    Vector2d AB(R(0,1) - R(0,0), R(1,1) - R(1,0));
    Vector2d BC(R(0,2) - R(0,1), R(1,2) - R(1,1));
    Vector2d CD(R(0,3) - R(0,2), R(1,3) - R(1,2));
    Vector2d DA(R(0,0) - R(0,3), R(1,0) - R(1,3));
    
    z[0] = cross2d(AB, BC);
    z[1] = cross2d(BC, CD);
    z[2] = cross2d(CD, DA);
    z[3] = cross2d(DA, AB);
    
    // Check if all cross products are nonzero and have the same sign
    bool all_positive = true, all_negative = true;
    for (int i = 0; i < 4; i++) {
        if (z[i] == 0.0) {
            throw runtime_error("Quadrilateral is not convex (collinear points)");
        }
        if (z[i] < 0.0) all_positive = false;
        if (z[i] > 0.0) all_negative = false;
    }
    if (!all_positive && !all_negative) {
        throw runtime_error("Quadrilateral is not convex");
    }
    
    // Mapping matrix M for triangles (Fortran data statement)
    // M = [2,4,1, 3,1,2, 4,2,3, 1,3,4] (12 elements, 4 triangles of 3 nodes each)
    // Convert to 0-based indices for internal use
    vector<vector<int>> triangle_nodes = {
        {1,3,0}, // triangle 1: B,D,A (nodes[1], nodes[3], nodes[0])
        {2,0,1}, // triangle 2: C,A,B (nodes[2], nodes[0], nodes[1])
        {3,1,2}, // triangle 3: D,B,C (nodes[3], nodes[1], nodes[2])
        {0,2,3}  // triangle 4: A,C,D (nodes[0], nodes[2], nodes[3])
    };
    
    // jnot = npivot + 2 if npivot <= 2, else npivot - 2 (1-based)
    int jnot = (npivot <= 2) ? npivot + 2 : npivot - 2;
    
    // Build full 12x12 stiffness matrix
    MatrixXd K_full = MatrixXd::Zero(12, 12);
    
    // Process each triangle except jnot (1-based index)
    for (size_t tri_idx = 0; tri_idx < triangle_nodes.size(); tri_idx++) {
        int tri_num_1based = static_cast<int>(tri_idx) + 1;
        if (tri_num_1based == jnot) continue;
        
        const vector<int>& tri = triangle_nodes[tri_idx];
        Vector3d p1 = nodes[tri[0]];
        Vector3d p2 = nodes[tri[1]];
        Vector3d p3 = nodes[tri[2]];
        
        MatrixXd K_super_u;
        try {
            ktrbsc(p1, p2, p3, youngs_modulus, nu, I, K_super_u);
        } catch (const exception& e) {
            // Re-throw with context
            ostringstream oss;
            oss << "Failed to compute triangle " << tri_num_1based << ": " << e.what();
            throw runtime_error(oss.str());
        }
        
        // Add contribution: each node has 3 DOF (w, theta_x, theta_y)
        // K_super_u is 9x9: blocks [0:2,0:2], [0:2,3:5], [0:2,6:8], [3:5,0:2], etc.
        // Map local node i (0,1,2) to global node tri[i]
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                int global_i = tri[i];
                int global_j = tri[j];
                int row_start = global_i * 3;
                int col_start = global_j * 3;
                
                // Extract 3x3 block: rows i*3 to i*3+2, cols j*3 to j*3+2
                MatrixXd block = K_super_u.block(i*3, j*3, 3, 3) / 2.0;
                K_full.block(row_start, col_start, 3, 3) += block;
            }
        }
    }
    
    return K_full;
}

int main() {
    try {
        // Test case: Quad A=(0,0,0), B=(1,0,0), C=(1,1,0), D=(0,1,0)
        vector<Vector3d> nodes = {
            Vector3d(0.0, 0.0, 0.0), // A
            Vector3d(1.0, 0.0, 0.0), // B
            Vector3d(1.0, 1.0, 0.0), // C
            Vector3d(0.0, 1.0, 0.0)  // D
        };
        
        double youngs_modulus = 200e9;      // Pa
        double nu = 0.3;       // Poisson's ratio
        double t = 0.01;       // m
        
        // Compute stiffness matrix
        MatrixXd K = kqdplt(nodes, youngs_modulus, nu, t);
        
        // Output as JSON-like array string
        cout << "{\"stiffness_matrix\":[" << endl;
        for (int i = 0; i < 12; i++) {
            cout << "  [";
            for (int j = 0; j < 12; j++) {
                cout << scientific << setprecision(15) << K(i, j);
                if (j < 11) cout << ",";
            }
            cout << "]";
            if (i < 11) cout << ",";
            cout << endl;
        }
        cout << "]}" << endl;
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}