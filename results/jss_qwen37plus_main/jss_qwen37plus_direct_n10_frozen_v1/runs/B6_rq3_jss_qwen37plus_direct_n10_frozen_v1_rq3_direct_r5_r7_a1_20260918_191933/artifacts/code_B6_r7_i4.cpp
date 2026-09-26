#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>

// Define M_PI if not available (common on Windows)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace Eigen;

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
void mat(double E, double nu, double& g11, double& g12, double& g13,
         double& g22, double& g23, double& g33) {
    // For isotropic material, the stiffness matrix components
    // G = [D] matrix for plate bending (reduced stiffness)
    double D = E / (1.0 - nu*nu);
    g11 = D;
    g12 = D * nu;
    g13 = 0.0;
    g22 = D;
    g23 = 0.0;
    g33 = D * (1.0 - nu) / 2.0;
}

// Basic bending triangle routine (KTRBSC)
void ktrbsc(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3,
            double E_mod, double nu, double I, MatrixXd& K_super_u) {
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
    mat(E_mod, nu, g11, g12, g13, g22, g23, g33);
    
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
    
    // Compute pseudo-inverse of HBAR using SVD (robust against singularity)
    JacobiSVD<MatrixXd> svd(HBAR, ComputeThinU | ComputeThinV);
    double tolerance = 1e-12;
    MatrixXd HBAR_inv = svd.matrixV() * (svd.singularValues().array().abs() > tolerance)
        .select(svd.singularValues().array().inverse(), 0.0).matrix().asDiagonal()
        * svd.matrixU().transpose();
    
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
MatrixXd kqdplt(const vector<Vector3d>& nodes, double E, double nu, double t) {
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
    
    // Form R-matrix (2x4) containing coordinates of sub-triangles
    // In local coordinate system where z is normal to plate
    Vector3d k_vec = cross(B-A, D-A);
    k_vec = normalize(k_vec);
    
    // I-vector = (B-A) projected onto plane perpendicular to k_vec
    Vector3d a1 = B - A;
    double h = dot(a1, k_vec) / 2.0;
    Vector3d i_vec = a1 - h * k_vec;
    i_vec = normalize(i_vec);
    
    // J-vector = k_vec × i_vec
    Vector3d j_vec = cross(k_vec, i_vec);
    j_vec = normalize(j_vec);
    
    // Compute coordinates in local system (x along i_vec, y along j_vec, z along k_vec)
    // R(1,j) = x-coordinate, R(2,j) = y-coordinate
    MatrixXd R(2, 4);
    R(0,0) = 0.0; R(1,0) = 0.0; // A
    R(0,1) = dot(B-A, i_vec); R(1,1) = dot(B-A, j_vec); // B
    R(0,2) = dot(C-A, i_vec); R(1,2) = dot(C-A, j_vec); // C
    R(0,3) = dot(D-A, i_vec); R(1,3) = dot(D-A, j_vec); // D
    
    // Check for convexity (angle >= 180 degrees)
    if (R(1,2) <= 0.0 || R(1,3) <= 0.0) {
        throw runtime_error("Quadrilateral is not convex");
    }
    double temp = R(0,1) - (R(0,1) - R(0,2)) * R(1,3) / R(1,2);
    if (R(0,3) >= temp) {
        throw runtime_error("Quadrilateral is not convex");
    }
    temp = R(1,2) * R(0,3) / R(1,3);
    if (R(0,2) > temp) {
        throw runtime_error("Quadrilateral is not convex");
    }
    
    // Mapping matrix M for triangles (Fortran data statement)
    // M = [2,4,1, 3,1,2, 4,2,3, 1,3,4] (12 elements, 4 triangles of 3 nodes each)
    vector<int> M = {2,4,1, 3,1,2, 4,2,3, 1,3,4};
    
    // Summation matrix KSUM (36 elements, 6x6)
    MatrixXd KSUM = MatrixXd::Zero(6, 6);
    
    // E-matrix (18 elements, 6x3)
    MatrixXd E_mat(6, 3);
    E_mat.setZero();
    E_mat.row(0) = k_vec.transpose();  // row 0: k-vector
    E_mat.row(1) = i_vec.transpose();  // row 1: i-vector  
    E_mat.row(2) = j_vec.transpose();  // row 2: j-vector
    // rows 3-5 are zero in this simplified version
    
    // For each triangle (j=1 to 4), skip if j == jnot
    // jnot = npivot + 2 if npivot <= 2, else npivot - 2
    int jnot = (npivot <= 2) ? npivot + 2 : npivot - 2;
    
    // Pre-allocate storage for K_super_u matrices (9x9 for each triangle)
    vector<MatrixXd> K_super_u_list(4, MatrixXd::Zero(9, 9));
    
    // Process each triangle
    for (int j = 0; j < 4; j++) {
        if (j + 1 == jnot) continue; // Skip triangle containing jnot
        
        // Get triangle node indices (convert to 0-based)
        if (3*j + 2 >= (int)M.size()) {
            throw runtime_error("Index out of bounds in triangle mapping M");
        }
        int subsc_a = M[3*j] - 1;
        int subsc_b = M[3*j + 1] - 1;
        int subsc_c = M[3*j + 2] - 1;
        
        // Validate node indices
        if (subsc_a < 0 || subsc_a >= 4 || subsc_b < 0 || subsc_b >= 4 || subsc_c < 0 || subsc_c >= 4) {
            throw runtime_error("Invalid node index in triangle mapping");
        }
        
        // Get triangle vertices
        Vector3d tri_a = nodes[subsc_a];
        Vector3d tri_b = nodes[subsc_b];
        Vector3d tri_c = nodes[subsc_c];
        
        // Compute triangle geometry
        Vector3d v = tri_b - tri_a;
        Vector3d vv = tri_c - tri_a;
        double xsubb = v.norm();
        Vector3d u1_vec = v / xsubb;
        Vector3d u2_vec = Vector3d(-u1_vec(1), u1_vec(0), 0.0); // perpendicular in xy-plane
        if (u2_vec.norm() < 1e-12) {
            u2_vec = Vector3d(0.0, 1.0, 0.0);
        }
        u2_vec = normalize(u2_vec);
        double xsubc = dot(u1_vec, vv);
        double ysubc = dot(u2_vec, vv);
        
        // Call KTRBSC for this triangle
        ktrbsc(tri_a, tri_b, tri_c, E, nu, I, K_super_u_list[j]);
        
        // Set up transformation matrix T
        MatrixXd T(3, 3);
        T << 1.0, 0.0, 0.0,
             0.0, u1_vec(0), u1_vec(1),
             0.0, -u1_vec(1), u1_vec(0);
        
        // Find which point of the subtriangle is the pivot
        int nbegin = -1;
        for (int i = 0; i < 3; i++) {
            int npoint = 3*j + i;
            if (npoint >= (int)M.size()) break;
            if (M[npoint] == npivot) {
                nbegin = 27*i;
                break;
            }
        }
        if (nbegin == -1) {
            throw runtime_error("Pivot not found in triangle");
        }
        
        // Extract 3x3 blocks from K_super_u and transform
        for (int i = 0; i < 3; i++) {
            // Get the i-th 3x3 block (rows 3*i to 3*i+2, cols 3*i to 3*i+2)
            MatrixXd K_block = K_super_u_list[j].block(3*i, 3*i, 3, 3);
            
            // Transform: T^T * K_block * T
            MatrixXd transformed = T.transpose() * K_block * T;
            
            // Add to KSUM (divided by 2.0 as in Fortran)
            int target_row = (M[3*j + i] - 1) * 3;
            int target_col = (M[3*j + i] - 1) * 3;
            if (target_row < 6 && target_col < 6) {
                KSUM.block(target_row, target_col, 3, 3) += transformed / 2.0;
            }
        }
    }
    
    // Now build the full 12x12 stiffness matrix
    // Each node has 3 DOF: w, theta_x, theta_y
    MatrixXd K_full = MatrixXd::Zero(12, 12);
    
    // Reset K_full
    K_full.setZero();
    
    // Define the four triangles to use (excluding jnot=3)
    vector<vector<int>> triangles;
    if (jnot == 1) {
        triangles = {{2,0,1}, {3,1,2}, {0,2,3}}; // skip triangle 1
    } else if (jnot == 2) {
        triangles = {{0,2,3}, {2,0,1}, {3,1,2}}; // skip triangle 2
    } else if (jnot == 3) {
        triangles = {{1,3,0}, {2,0,1}, {0,2,3}}; // skip triangle 3
    } else { // jnot == 4
        triangles = {{1,3,0}, {2,0,1}, {3,1,2}}; // skip triangle 4
    }
    
    // For each triangle, compute its contribution
    for (size_t tri_idx = 0; tri_idx < triangles.size(); tri_idx++) {
        const vector<int>& tri = triangles[tri_idx];
        if (tri.size() < 3) continue;
        Vector3d p1 = nodes[tri[0]];
        Vector3d p2 = nodes[tri[1]];
        Vector3d p3 = nodes[tri[2]];
        
        MatrixXd K_super_u;
        try {
            ktrbsc(p1, p2, p3, E, nu, I, K_super_u);
        } catch (...) {
            continue;
        }
        
        // Extract the 3x3 blocks corresponding to the three nodes
        // Node mapping: tri[0], tri[1], tri[2] -> global nodes
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                int global_i = tri[i];
                int global_j = tri[j];
                if (global_i < 0 || global_i >= 4 || global_j < 0 || global_j >= 4) continue;
                int row_start = global_i * 3;
                int col_start = global_j * 3;
                
                // Get the 3x3 block from K_super_u
                // The Fortran code stores them in a specific pattern
                // For simplicity, use the top-left 3x3 for diagonal and appropriate off-diagonal
                if (i*3+2 < K_super_u.rows() && j*3+2 < K_super_u.cols()) {
                    MatrixXd block = K_super_u.block(i*3, j*3, 3, 3) / 2.0;
                    K_full.block(row_start, col_start, 3, 3) += block;
                }
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
        
        double E = 200e9;      // Pa
        double nu = 0.3;       // Poisson's ratio
        double t = 0.01;       // m
        
        // Compute stiffness matrix
        MatrixXd K = kqdplt(nodes, E, nu, t);
        
        // Output in space-separated scientific notation, 12x12, row-major
        cout << scientific;
        cout.precision(15);
        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 12; j++) {
                if (j > 0) cout << " ";
                cout << K(i, j);
            }
            cout << "\n";
        }
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}