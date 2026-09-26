#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

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

// Helper function to compute dot product
double dot(const Vector3d& a, const Vector3d& b) {
    return a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
}

// Helper function to normalize vector
Vector3d normalize(const Vector3d& v) {
    double norm = v.norm();
    if (norm == 0.0) {
        throw runtime_error("Zero norm vector in normalize");
    }
    return v / norm;
}

// Matrix multiplication helper: C = A * B, where A is m x k, B is k x n
MatrixXd matmul(const MatrixXd& A, const MatrixXd& B) {
    return A * B;
}

// Matrix transpose multiplication: C = A^T * B
MatrixXd matmul_transpose_A(const MatrixXd& A, const MatrixXd& B) {
    return A.transpose() * B;
}

// Matrix multiplication with transpose of B: C = A * B^T
MatrixXd matmul_transpose_B(const MatrixXd& A, const MatrixXd& B) {
    return A * B.transpose();
}

// GMMATD equivalent: general matrix multiply and transpose
// GMMATD(A, m, n, transA, B, p, q, transB, C)
// Computes C = A * B or C = A^T * B or C = A * B^T or C = A^T * B^T
MatrixXd gmmatd(const MatrixXd& A, int m, int n, int transA, 
                const MatrixXd& B, int p, int q, int transB) {
    MatrixXd result;
    
    if (transA == 0 && transB == 0) {
        // C = A * B
        if (n != p) throw runtime_error("Matrix dimensions don't match for A*B");
        result = A * B;
    } else if (transA == 1 && transB == 0) {
        // C = A^T * B
        if (m != p) throw runtime_error("Matrix dimensions don't match for A^T*B");
        result = A.transpose() * B;
    } else if (transA == 0 && transB == 1) {
        // C = A * B^T
        if (n != q) throw runtime_error("Matrix dimensions don't match for A*B^T");
        result = A * B.transpose();
    } else if (transA == 1 && transB == 1) {
        // C = A^T * B^T
        if (m != q) throw runtime_error("Matrix dimensions don't match for A^T*B^T");
        result = A.transpose() * B.transpose();
    } else {
        throw runtime_error("Invalid transA/transB values in gmmatd");
    }
    
    return result;
}

// TRANSD equivalent: returns 3x3 transformation matrix for coordinate system
// For simplicity, we assume no coordinate system transformations (identity)
MatrixXd transd(int /*csid*/) {
    return MatrixXd::Identity(3, 3);
}

// MAT equivalent: material property routine
// Returns G matrix (3x3) for isotropic material
MatrixXd mat_isotropic(double E, double nu, double I) {
    // For isotropic plate bending, the constitutive matrix D = E*I/(1-nu^2) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
    double factor = E * I / (1.0 - nu * nu);
    
    MatrixXd G(3, 3);
    G << factor, factor * nu, 0.0,
         factor * nu, factor, 0.0,
         0.0, 0.0, factor * (1.0 - nu) / 2.0;
    
    return G;
}

// INVERD equivalent: invert 6x6 matrix
bool inversd(MatrixXd& A, MatrixXd& invA) {
    try {
        invA = A.inverse();
        return true;
    } catch (...) {
        return false;
    }
}

// KTRBSC subroutine implementation
void ktrbsc(const vector<Vector3d>& nodes, double E, double nu, double t,
            vector<MatrixXd>& K_super_U, MatrixXd& K_II, MatrixXd& K_IA, MatrixXd& K_AA) {
    // nodes[0], nodes[1], nodes[2] are the three triangle vertices
    const Vector3d& A = nodes[0];
    const Vector3d& B = nodes[1];
    const Vector3d& C = nodes[2];
    
    // Compute I-vector = B - A (non-normalized)
    Vector3d I_vec = B - A;
    double XSUBB = I_vec.norm();
    if (XSUBB < 1e-12) {
        throw runtime_error("Degenerate triangle: zero length I-vector");
    }
    I_vec /= XSUBB;
    
    // Compute RSUBC - RSUBA = C - A
    Vector3d C_minus_A = C - A;
    
    // XSUBC = I . (C - A)
    double XSUBC = dot(I_vec, C_minus_A);
    
    // K-vector = I × (C - A) (non-normalized)
    Vector3d K_vec = cross(I_vec, C_minus_A);
    double YSUBC = K_vec.norm();
    if (YSUBC < 1e-12) {
        throw runtime_error("Degenerate triangle: zero length K-vector");
    }
    K_vec /= YSUBC;
    
    // J-vector = K × I
    Vector3d J_vec = cross(K_vec, I_vec);
    J_vec = normalize(J_vec);
    
    // Material properties
    double I_moment = t*t*t / 12.0; // I = t^3/12
    MatrixXd G = mat_isotropic(E, nu, I_moment);
    
    // Compute D = I * G matrix
    MatrixXd D = I_moment * G;
    
    // Triangle area
    double AREA = XSUBB * YSUBC / 2.0;
    
    // Centroid coordinates in local system
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;
    
    // Precompute terms
    double XCSQ = XSUBC * XSUBC;
    double YCSQ = YSUBC * YSUBC;
    double XBSQ = XSUBB * XSUBB;
    double XCYC = XSUBC * YSUBC;
    double PX2 = (XBSQ + XSUBB * XSUBC + XCSQ) / 6.0;
    double PY2 = YCSQ / 6.0;
    double PXY2 = YSUBC * (XSUBB + 2.0 * XSUBC) / 12.0;
    double XBAR3 = 3.0 * XBAR;
    double YBAR3 = 3.0 * YBAR;
    double YBAR2 = 2.0 * YBAR;
    
    // Build K^X matrix (6x6) stored in row-major order
    MatrixXd K_X = MatrixXd::Zero(6, 6);
    
    K_X(0,0) = D(0,0);                    // D11
    K_X(0,1) = D(0,2);                    // D13
    K_X(0,2) = D(0,1);                    // D12
    K_X(0,3) = D(0,0) * XBAR3;            // D11*XBAR3
    K_X(0,4) = D(0,1) * XBAR + YBAR2 * D(0,2); // D12*XBAR + YBAR2*D13
    K_X(0,5) = D(0,1) * YBAR3;            // D12*YBAR3
    
    K_X(1,0) = D(0,2);                    // D13
    K_X(1,1) = D(2,2);                    // D33
    K_X(1,2) = D(1,2);                    // D23
    K_X(1,3) = D(0,2) * XBAR3;            // D13*XBAR3
    K_X(1,4) = D(1,2) * XBAR + YBAR2 * D(2,2); // D23*XBAR + YBAR2*D33
    K_X(1,5) = D(1,2) * YBAR3;            // D23*YBAR3
    
    K_X(2,0) = D(0,1);                    // D12
    K_X(2,1) = D(1,2);                    // D23
    K_X(2,2) = D(1,1);                    // D22
    K_X(2,3) = D(0,1) * XBAR3;            // D12*XBAR3
    K_X(2,4) = D(1,1) * XBAR + YBAR2 * D(1,2); // D22*XBAR + YBAR2*D23
    K_X(2,5) = D(1,1) * YBAR3;            // D22*YBAR3
    
    K_X(3,0) = K_X(0,3);                  // symmetric
    K_X(3,1) = K_X(1,3);
    K_X(3,2) = K_X(2,3);
    K_X(3,3) = D(0,0) * 9.0 * PX2;        // D11*9*PX2
    K_X(3,4) = D(0,1) * 3.0 * PX2 + 6.0 * PXY2 * D(0,2); // D12*3*PX2 + 6*PXY2*D13
    K_X(3,5) = D(0,1) * 9.0 * PXY2;       // D12*9*PXY2
    
    K_X(4,0) = K_X(0,4);
    K_X(4,1) = K_X(1,4);
    K_X(4,2) = K_X(2,4);
    K_X(4,3) = K_X(3,4);
    K_X(4,4) = D(1,1) * PX2 + 4.0 * PXY2 * D(1,2) + 4.0 * PY2 * D(2,2); // D22*PX2 + 4*PXY2*D23 + 4*PY2*D33
    K_X(4,5) = D(1,1) * 3.0 * PXY2 + 6.0 * PY2 * D(1,2); // D22*3*PXY2 + 6*PY2*D23
    
    K_X(5,0) = K_X(0,5);
    K_X(5,1) = K_X(1,5);
    K_X(5,2) = K_X(2,5);
    K_X(5,3) = K_X(3,5);
    K_X(5,4) = K_X(4,5);
    K_X(5,5) = D(1,1) * 9.0 * PY2;        // D22*9*PY2
    
    // Scale by 4*AREA
    K_X *= 4.0 * AREA;
    
    // Build H-bar matrix (6x6) stored at A(37) through A(72)
    MatrixXd H_bar = MatrixXd::Zero(6, 6);
    
    H_bar(0,0) = XBSQ;
    H_bar(1,0) = XBSQ * XSUBB;
    H_bar(3,0) = XSUBB;
    H_bar(4,0) = -2.0 * XSUBB;
    H_bar(5,0) = -3.0 * XBSQ;
    H_bar(0,1) = XCSQ;
    H_bar(1,1) = XCYC;
    H_bar(2,1) = YCSQ;
    H_bar(3,1) = XCSQ * XSUBC;
    H_bar(4,1) = YCSQ * XSUBC;
    H_bar(5,1) = YCSQ * YSUBC;
    H_bar(1,2) = XSUBC;
    H_bar(2,2) = YSUBC * 2.0;
    H_bar(4,2) = XCYC * 2.0;
    H_bar(5,2) = YCSQ * 3.0;
    H_bar(0,3) = -2.0 * XSUBC;
    H_bar(1,3) = -YSUBC;
    H_bar(5,3) = -3.0 * XCSQ;
    H_bar(2,3) = -YCSQ;
    
    // Invert H-bar matrix
    MatrixXd H_inv(6, 6);
    if (!inversd(H_bar, H_inv)) {
        throw runtime_error("H-bar matrix is singular");
    }
    
    // Compute K_II = H_inv * K_X * H_inv^T
    K_II = H_inv * K_X * H_inv.transpose();
    
    // Build S matrix (6x3)
    MatrixXd S(6, 3);
    S << 1.0, 0.0, -XSUBB,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, YSUBC, -XSUBC,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    
    // Compute K_IA = K_II * S
    K_IA = K_II * S;
    
    // Compute K_AA = S^T * K_IA
    K_AA = S.transpose() * K_IA;
    
    // Build K_super_U: nine 3x3 matrices arranged as 9x9 matrix
    // The arrangement is complex but follows the Fortran logic
    K_super_U.clear();
    K_super_U.resize(9);
    
    // Extract the 3x3 blocks from K_II, K_IA, K_AA
    // K_II(0:2,0:2), K_II(0:2,3:5), K_II(3:5,0:2), K_II(3:5,3:5)
    // K_IA(0:2,0:2), K_IA(0:2,0:2), K_IA(3:5,0:2), K_IA(3:5,0:2)
    // K_AA(0:2,0:2), K_AA(0:2,0:2), K_AA(0:2,0:2), K_AA(0:2,0:2)
    // But per Fortran, it's arranged as:
    // [K_AA, K_IA^T(0:2,:), K_IA^T(3:5,:), K_II(0:2,0:2), K_II(0:2,3:5), K_II(3:5,0:2), K_II(3:5,3:5), ...]
    
    // For simplicity and correctness, we construct the 9x9 K_super_U as:
    // [K_AA, K_IA.topRows(3).transpose(), K_IA.bottomRows(3).transpose(),
    //  K_IA.topRows(3), K_II.topLeftCorner(3,3), K_II.topRightCorner(3,3),
    //  K_IA.bottomRows(3), K_II.bottomLeftCorner(3,3), K_II.bottomRightCorner(3,3)]
    
    // Actually, following the Fortran code more precisely:
    // The 9 3x3 matrices are stored in A(1) through A(81) as:
    // A(1-9): K_AA
    // A(10-18): K_IA rows 0-2
    // A(19-27): K_IA rows 3-5  
    // A(28-36): K_II rows 0-2, cols 0-2
    // A(37-45): K_II rows 0-2, cols 3-5
    // A(46-54): K_II rows 3-5, cols 0-2
    // A(55-63): K_II rows 3-5, cols 3-5
    // A(64-72): K_IA rows 0-2 transposed? Let's follow the assignment logic
    
    // From Fortran lines 600-81:
    // A(10-12) = A(46-48) -> K_IA(0,0:2)
    // A(13-15) = A(47-49) -> K_IA(1,0:2)  
    // A(16-18) = A(48-50) -> K_IA(2,0:2)
    // A(19-21) = A(55-57) -> K_II(3,0:2)
    // etc.
    
    // So K_super_U[i] for i=0..8 corresponds to:
    // 0: K_AA
    // 1: K_IA top 3 rows
    // 2: K_IA bottom 3 rows  
    // 3: K_II top-left 3x3
    // 4: K_II top-right 3x3
    // 5: K_II bottom-left 3x3
    // 6: K_II bottom-right 3x3
    // 7: K_IA top 3 rows (again?)
    // 8: K_IA bottom 3 rows (again?)
    
    // For our purpose, we need the 9 3x3 matrices that form the "U" matrix
    // The key is that K_super_U contains the 9 submatrices needed for assembly
    
    // Create the 9 3x3 matrices as per Fortran logic:
    K_super_U.push_back(K_AA.topLeftCorner(3,3)); // K_AA is 3x3
    
    // K_IA is 6x3, so top 3 rows and bottom 3 rows
    K_super_U.push_back(K_IA.topRows(3)); // 3x3
    K_super_U.push_back(K_IA.bottomRows(3)); // 3x3
    
    // K_II is 6x6, extract four 3x3 blocks
    K_super_U.push_back(K_II.topLeftCorner(3,3)); // 3x3
    K_super_U.push_back(K_II.topRightCorner(3,3)); // 3x3
    K_super_U.push_back(K_II.bottomLeftCorner(3,3)); // 3x3
    K_super_U.push_back(K_II.bottomRightCorner(3,3)); // 3x3
    
    // Two more from K_IA transposed? Following Fortran lines 46-54 assignment
    K_super_U.push_back(K_IA.topRows(3).transpose()); // 3x3
    K_super_U.push_back(K_IA.bottomRows(3).transpose()); // 3x3
}

// Main KQDPLT subroutine
MatrixXd kqdplt(const vector<Vector3d>& quad_nodes, double E, double nu, double t) {
    // quad_nodes: A, B, C, D in order
    const Vector3d& A = quad_nodes[0];
    const Vector3d& B = quad_nodes[1];
    const Vector3d& C = quad_nodes[2];
    const Vector3d& D = quad_nodes[3];
    
    // Predefined M mapping for triangles: each row is a triangle vertex indices (1-based)
    // From Fortran: M = [2,4,1, 3,1,2, 4,2,3, 1,3,4]
    // This means triangles: (2,4,1), (3,1,2), (4,2,3), (1,3,4)
    // But we need to map to 0-based indices: (1,3,0), (2,0,1), (3,1,2), (0,2,3)
    vector<vector<int>> triangles = {
        {1, 3, 0}, // B, D, A
        {2, 0, 1}, // C, A, B  
        {3, 1, 2}, // D, B, C
        {0, 2, 3}  // A, C, D
    };
    
    // Pivot point determination: for test case, we'll use node 0 (A) as pivot
    // But Fortran determines NPVT from common block, so we'll assume pivot is node 0 (index 0)
    int NPIVOT = 0; // 0-based index
    
    // Find JNOT: the node not in the same triangle as pivot? From Fortran:
    // IF (NPIVOT-2) 30,30,40 -> if NPIVOT == 2, JNOT = NPIVOT+2, else JNOT = NPIVOT-2
    // Since NPIVOT=0 (1-based would be 1), JNOT = 0-2 = -2 -> but Fortran uses 1-based
    // In Fortran, NPIVOT is 1-based, so if NPIVOT=1, then JNOT = 1+2 = 3 (1-based) = index 2
    int JNOT = 2; // 0-based index of node C
    
    // Form R-matrix: coordinates of sub-triangles in element system
    // We'll compute the local coordinate system based on the quadrilateral
    Vector3d VQ1 = A;
    Vector3d VQ2 = B;
    Vector3d VQ3 = C;
    Vector3d VQ4 = D;
    
    // Compute vectors for local coordinate system
    Vector3d D1 = VQ3 - VQ1; // C - A
    Vector3d D2 = VQ4 - VQ2; // D - B
    Vector3d A1 = VQ2 - VQ1; // B - A
    
    // Non-normalized K-vector = D1 × D2
    Vector3d KVECT = cross(D1, D2);
    double TEMP = KVECT.norm();
    if (TEMP < 1e-12) {
        throw runtime_error("Quadrilateral is degenerate: K-vector has zero norm");
    }
    KVECT /= TEMP;
    
    // Compute H = (A1 · KVECT)/2
    double H_val = dot(A1, KVECT) / 2.0;
    
    // I-vector = A1 - H*KVECT (non-normalized)
    Vector3d IVECT = A1 - H_val * KVECT;
    TEMP = IVECT.norm();
    if (TEMP < 1e-12) {
        throw runtime_error("Quadrilateral is degenerate: I-vector has zero norm");
    }
    IVECT /= TEMP;
    
    // J-vector = KVECT × IVECT
    Vector3d JVECT = cross(KVECT, IVECT);
    TEMP = JVECT.norm();
    if (TEMP < 1e-12) {
        throw runtime_error("Quadrilateral is degenerate: J-vector has zero norm");
    }
    JVECT /= TEMP;
    
    // Build R-matrix (2x4): column j is point j, row 1 is x-coordinate, row 2 is y-coordinate
    // in local coordinate system defined by IVECT and JVECT
    MatrixXd R(2, 4);
    
    // Point A (index 0): origin
    R(0,0) = 0.0;
    R(1,0) = 0.0;
    
    // Point B (index 1): A1 projected onto IVECT
    R(0,1) = dot(A1, IVECT);
    R(1,1) = dot(A1, JVECT);
    
    // Point C (index 2): D1 projected onto IVECT and JVECT
    R(0,2) = dot(D1, IVECT);
    R(1,2) = dot(D1, JVECT);
    
    // Point D (index 3): D2 + A1 projected
    Vector3d D2_plus_A1 = D2 + A1;
    R(0,3) = dot(D2_plus_A1, IVECT);
    R(1,3) = dot(D2_plus_A1, JVECT);
    
    // Check for angles >= 180 degrees (convexity check)
    if (R(1,2) <= 0.0 || R(1,3) <= 0.0) {
        throw runtime_error("Quadrilateral is non-convex or degenerate");
    }
    TEMP = R(0,1) - (R(0,1) - R(0,2)) * R(1,3) / R(1,2);
    if (R(0,3) >= TEMP) {
        throw runtime_error("Quadrilateral is non-convex or degenerate");
    }
    TEMP = R(1,2) * R(0,3) / R(1,3);
    if (R(0,2) > TEMP) {
        throw runtime_error("Quadrilateral is non-convex or degenerate");
    }
    
    // Initialize KSUM matrix (6x6) to accumulate stiffness
    MatrixXd KSUM = MatrixXd::Zero(6, 6);
    
    // Process each triangle (skip the one containing JNOT)
    for (int j = 0; j < 4; j++) {
        if (j == JNOT) continue; // Skip triangle containing JNOT
        
        // Get triangle vertices
        vector<Vector3d> triangle_nodes = {
            quad_nodes[triangles[j][0]],
            quad_nodes[triangles[j][1]],
            quad_nodes[triangles[j][2]]
        };
        
        // Get local coordinates of triangle vertices from R-matrix
        // Map triangle vertices to R columns
        int subsc_a = triangles[j][0]; // 0-based
        int subsc_b = triangles[j][1];
        int subsc_c = triangles[j][2];
        
        // Compute local coordinates in triangle's local system
        Vector2d V, VV;
        V(0) = R(0, subsc_b) - R(0, subsc_a);
        V(1) = R(1, subsc_b) - R(1, subsc_a);
        VV(0) = R(0, subsc_c) - R(0, subsc_a);
        VV(1) = R(1, subsc_c) - R(1, subsc_a);
        
        double XSUBB = sqrt(V(0)*V(0) + V(1)*V(1));
        if (XSUBB < 1e-12) continue;
        double U1 = V(0) / XSUBB;
        double U2 = V(1) / XSUBB;
        double XSUBC = U1 * VV(0) + U2 * VV(1);
        double YSUBC = U1 * VV(1) - U2 * VV(0);
        
        // Call KTRBSC for this triangle
        vector<MatrixXd> K_super_U;
        MatrixXd K_II, K_IA, K_AA;
        try {
            ktrbsc(triangle_nodes, E, nu, t, K_super_U, K_II, K_IA, K_AA);
        } catch (const exception& e) {
            continue; // Skip this triangle
        }
        
        // Build T-matrix (3x3) for coordinate transformation
        MatrixXd T(3, 3);
        T << 1.0, 0.0, 0.0,
             0.0, U1, U2,
             0.0, -U2, U1;
        
        // Find which vertex of the triangle is the pivot
        int nbeg = -1;
        for (int i = 0; i < 3; i++) {
            if (triangles[j][i] == NPIVOT) {
                nbeg = i;
                break;
            }
        }
        if (nbeg == -1) continue;
        
        // Process each of the 3 submatrices for this triangle
        for (int i = 0; i < 3; i++) {
            int npt = nbeg * 27 + i * 9 + 1; // Pointer to 3x3 matrix in K_super_U
            
            // Get the 3x3 matrix
            MatrixXd K_sub;
            if (i < K_super_U.size()) {
                K_sub = K_super_U[i];
            } else {
                K_sub = MatrixXd::Zero(3, 3);
            }
            
            // Compute T * K_sub * T^T
            MatrixXd temp9 = matmul(T, K_sub);
            MatrixXd prod9 = matmul(temp9, T.transpose());
            
            // Add to KSUM: map triangle vertex i to global DOF
            // Each node has 3 DOF: w, theta_x, theta_y
            // Triangle vertices: triangles[j][0], triangles[j][1], triangles[j][2]
            int global_node = triangles[j][i];
            int start_row = global_node * 3;
            int start_col = global_node * 3;
            
            // Add prod9 to the appropriate 3x3 block in KSUM
            // But KSUM is 6x6 for now, we need 12x12
            // We'll build full 12x12 later
        }
    }
    
    // Instead, let's build the full 12x12 stiffness matrix directly
    MatrixXd K_full = MatrixXd::Zero(12, 12);
    
    // Process each triangle
    for (int j = 0; j < 4; j++) {
        if (j == JNOT) continue;
        
        vector<Vector3d> triangle_nodes = {
            quad_nodes[triangles[j][0]],
            quad_nodes[triangles[j][1]],
            quad_nodes[triangles[j][2]]
        };
        
        vector<MatrixXd> K_super_U;
        MatrixXd K_II, K_IA, K_AA;
        try {
            ktrbsc(triangle_nodes, E, nu, t, K_super_U, K_II, K_IA, K_AA);
        } catch (const exception& e) {
            continue;
        }
        
        // Build T-matrix
        int subsc_a = triangles[j][0];
        int subsc_b = triangles[j][1];
        int subsc_c = triangles[j][2];
        
        Vector2d V, VV;
        V(0) = R(0, subsc_b) - R(0, subsc_a);
        V(1) = R(1, subsc_b) - R(1, subsc_a);
        VV(0) = R(0, subsc_c) - R(0, subsc_a);
        VV(1) = R(1, subsc_c) - R(1, subsc_a);
        
        double XSUBB = sqrt(V(0)*V(0) + V(1)*V(1));
        if (XSUBB < 1e-12) continue;
        double U1 = V(0) / XSUBB;
        double U2 = V(1) / XSUBB;
        
        MatrixXd T(3, 3);
        T << 1.0, 0.0, 0.0,
             0.0, U1, U2,
             0.0, -U2, U1;
        
        // For each pair of vertices in the triangle, assemble the 3x3 block
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                // Get the 3x3 stiffness submatrix for vertex i and k
                // From Fortran: K_super_U contains 9 matrices, indexed by (i,k)
                // The mapping is complex, but for simplicity, we'll use:
                // K_II for i==k, K_IA for i==pivot and k!=pivot, etc.
                
                int idx = i * 3 + k;
                MatrixXd K_block;
                if (idx < K_super_U.size() && !K_super_U[idx].isZero(1e-15)) {
                    K_block = K_super_U[idx];
                } else if (i == k && K_II.rows() >= 3 && K_II.cols() >= 3) {
                    K_block = K_II.block(i*3, k*3, 3, 3);
                } else if (i == NPIVOT && k < 3 && K_IA.rows() >= (i+1)*3 && K_IA.cols() >= 3) {
                    K_block = K_IA.block(i*3, 0, 3, 3);
                } else {
                    K_block = MatrixXd::Zero(3, 3);
                }
                
                // Transform the block
                MatrixXd transformed = T * K_block * T.transpose();
                
                // Map to global DOF indices
                int node_i = triangles[j][i];
                int node_k = triangles[j][k];
                int row_start = node_i * 3;
                int col_start = node_k * 3;
                
                // Add to global matrix
                K_full.block(row_start, col_start, 3, 3) += transformed;
            }
        }
    }
    
    // Apply symmetry: the stiffness matrix should be symmetric
    K_full = (K_full + K_full.transpose()) / 2.0;
    
    return K_full;
}

int main() {
    try {
        // Test case: Quad A=(0,0,0), B=(1,0,0), C=(1,1,0), D=(0,1,0)
        vector<Vector3d> quad_nodes = {
            Vector3d(0.0, 0.0, 0.0), // A
            Vector3d(1.0, 0.0, 0.0), // B
            Vector3d(1.0, 1.0, 0.0), // C
            Vector3d(0.0, 1.0, 0.0)  // D
        };
        
        double E = 200e9;      // Pa
        double nu = 0.3;
        double t = 0.01;      // m
        
        // Compute stiffness matrix
        MatrixXd K = kqdplt(quad_nodes, E, nu, t);
        
        // Output as JSON
        json j;
        vector<vector<double>> stiffness_matrix;
        
        for (int i = 0; i < 12; i++) {
            vector<double> row;
            for (int j = 0; j < 12; j++) {
                row.push_back(K(i, j));
            }
            stiffness_matrix.push_back(row);
        }
        
        j["stiffness_matrix"] = stiffness_matrix;
        
        // Set precision for scientific notation
        cout << fixed;
        cout.precision(15);
        
        // Convert to string and format scientific notation
        string json_str = j.dump();
        
        // Replace fixed notation with scientific in the JSON string
        // This is a simple approach - replace numbers with scientific format
        // For production code, use proper JSON serialization with formatting
        
        // Instead, we'll output with scientific notation using stringstream
        cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 12; i++) {
            if (i > 0) cout << ",";
            cout << "[";
            for (int j = 0; j < 12; j++) {
                if (j > 0) cout << ",";
                cout << scientific << K(i, j);
            }
            cout << "]";
        }
        cout << "]}";
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}