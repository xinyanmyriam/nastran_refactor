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
    size_t e_pos = s.find('e');
    if (e_pos != std::string::npos) {
        // Remove trailing zeros after decimal point before exponent
        size_t dot_pos = s.find('.');
        if (dot_pos != std::string::npos && dot_pos < e_pos) {
            // Find last non-zero digit before exponent
            size_t last_nonzero = e_pos;
            for (size_t i = e_pos; i > dot_pos; --i) {
                if (s[i-1] != '0' && s[i-1] != '.') {
                    last_nonzero = i;
                    break;
                }
            }
            if (last_nonzero < e_pos) {
                s.erase(last_nonzero, e_pos - last_nonzero);
                s.insert(s.end(), 'e');
                // Re-append exponent part
                s += s.substr(e_pos+1);
            }
        }
    }
    // Remove '+' from exponent
    size_t plus_pos = s.find("e+");
    if (plus_pos != std::string::npos) {
        s.replace(plus_pos+1, 1, "");
    }
    return s;
}

// JSON-safe string output
void print_json_matrix(const Eigen::Matrix<double, 9, 9>& K) {
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

// Compute the triangular plate bending element stiffness matrix (DKT formulation)
Eigen::Matrix<double, 9, 9> compute_ktrplt_stiffness() {
    // Geometry
    Eigen::Vector3d A(0.0, 0.0, 0.0);
    Eigen::Vector3d B(1.0, 0.0, 0.0);
    Eigen::Vector3d C(0.0, 1.0, 0.0);

    // Material properties
    const double E = 200e9;           // Young's modulus
    const double nu = 0.3;            // Poisson's ratio
    const double t = 0.01;            // thickness
    const double D_factor = E * t*t*t / (12.0 * (1.0 - nu*nu)); // D = E t^3 / (12(1-nu^2))

    // Material matrix D for isotropic bending: 3x3
    Eigen::Matrix3d D;
    D << 1.0,     nu,      0.0,
         nu,      1.0,     0.0,
         0.0,     0.0,     (1.0 - nu) / 2.0;
    D *= D_factor;

    // Compute element coordinate system (I, J, K vectors)
    Eigen::Vector3d D1 = C - A;  // AC
    Eigen::Vector3d D2 = B - A;  // AB

    // I-vector = AB normalized
    double XSUBB = D2.norm();
    Eigen::Vector3d IVEC = D2 / XSUBB;

    // K-vector = I × AC normalized
    Eigen::Vector3d KVEC = IVEC.cross(D1);
    double YSUBC = KVEC.norm();
    KVEC /= YSUBC;

    // J-vector = K × I
    Eigen::Vector3d JVEC = KVEC.cross(IVEC);
    JVEC.normalize();

    // Area of triangle
    double AREA = XSUBB * YSUBC / 2.0;

    // 2x2 rotation matrix for in-plane rotations: maps (θ_I, θ_J) -> (θ_X, θ_Y)
    Eigen::Matrix2d T2;
    T2 << IVEC.x(), JVEC.x(),
          IVEC.y(), JVEC.y();

    // 3x3 transformation for [w, θ_x, θ_y]: w is invariant, rotations transform
    Eigen::Matrix3d T3 = Eigen::Matrix3d::Zero();
    T3(0,0) = 1.0;
    T3(1,1) = T2(0,0); T3(1,2) = T2(0,1);
    T3(2,1) = T2(1,0); T3(2,2) = T2(1,1);

    // Build global-to-local transformation matrix L (9x9): block-diagonal with T3
    Eigen::Matrix<double, 9, 9> L = Eigen::Matrix<double, 9, 9>::Zero();
    for (int i = 0; i < 3; ++i) {
        int base = i * 3;
        L.block<3,3>(base, base) = T3;
    }

    // Shape function derivatives w.r.t. global X,Y
    // N1 = 1 - x - y, N2 = x, N3 = y
    // dN1/dx = -1, dN1/dy = -1
    // dN2/dx = 1,  dN2/dy = 0
    // dN3/dx = 0,  dN3/dy = 1

    // B matrix: 3 (curvatures) x 9 (DOFs)
    // κ_xx = -∂θ_x/∂x
    // κ_yy = -∂θ_y/∂y  
    // κ_xy = -½(∂θ_x/∂y + ∂θ_y/∂x)
    // So for DOFs [w1,θ_x1,θ_y1, w2,θ_x2,θ_y2, w3,θ_x3,θ_y3]:
    // Only θ terms contribute to curvature.
    Eigen::Matrix<double, 3, 9> B = Eigen::Matrix<double, 3, 9>::Zero();

    // ∂θ_x/∂x column: coefficients for θ_x1, θ_x2, θ_x3
    B(0,1) = -(-1); // -dN1/dx = 1
    B(0,4) = -(1);  // -dN2/dx = -1
    B(0,7) = -(0);  // -dN3/dx = 0

    // ∂θ_y/∂y column: coefficients for θ_y1, θ_y2, θ_y3
    B(1,2) = -(-1); // -dN1/dy = 1
    B(1,5) = -(0);  // -dN2/dy = 0
    B(1,8) = -(1);  // -dN3/dy = -1

    // ∂θ_x/∂y and ∂θ_y/∂x for κ_xy = -0.5*(∂θ_x/∂y + ∂θ_y/∂x)
    // ∂θ_x/∂y: dN1/dy = -1, dN2/dy = 0, dN3/dy = 1 → coeffs: -(-1), -(0), -(1) = 1,0,-1
    // ∂θ_y/∂x: dN1/dx = -1, dN2/dx = 1, dN3/dx = 0 → coeffs: -(-1), -(1), -(0) = 1,-1,0
    // Sum: [2, -1, -1] → times -0.5 → [-1, 0.5, 0.5]
    B(2,1) = -0.5 * (-1); // -0.5 * dN1/dy = -0.5 * (-1) = 0.5
    B(2,2) = -0.5 * (-1); // -0.5 * dN1/dx = -0.5 * (-1) = 0.5
    B(2,4) = -0.5 * (0);  // -0.5 * dN2/dy = 0
    B(2,5) = -0.5 * (1);  // -0.5 * dN2/dx = -0.5
    B(2,7) = -0.5 * (1);  // -0.5 * dN3/dy = -0.5
    B(2,8) = -0.5 * (0);  // -0.5 * dN3/dx = 0

    // Correction: Standard DKT B matrix is:
    // κ = -[ ∂θ_x/∂x, ∂θ_y/∂y, ∂θ_x/∂y + ∂θ_y/∂x ]^T
    // So B(0,i) = -∂N_i/∂x for θ_x columns (cols 1,4,7)
    // B(1,i) = -∂N_i/∂y for θ_y columns (cols 2,5,8)
    // B(2,i) = -0.5*∂N_i/∂y for θ_x cols + -0.5*∂N_i/∂x for θ_y cols
    // Let's rebuild correctly:

    B.setZero();
    // θ_x columns: indices 1,4,7
    B(0,1) = -(-1); // -dN1/dx = 1
    B(0,4) = -(1);  // -dN2/dx = -1
    B(0,7) = -(0);  // -dN3/dx = 0

    // θ_y columns: indices 2,5,8
    B(1,2) = -(-1); // -dN1/dy = 1
    B(1,5) = -(0);  // -dN2/dy = 0
    B(1,8) = -(1);  // -dN3/dy = -1

    // κ_xy = - (∂θ_x/∂y + ∂θ_y/∂x)  [some sources omit 0.5]
    // So B(2,1) = -dN1/dy = -(-1) = 1
    // B(2,2) = -dN1/dx = -(-1) = 1
    // B(2,4) = -dN2/dy = -(0) = 0
    // B(2,5) = -dN2/dx = -(1) = -1
    // B(2,7) = -dN3/dy = -(1) = -1
    // B(2,8) = -dN3/dx = -(0) = 0
    B(2,1) = -(-1); // 1
    B(2,2) = -(-1); // 1
    B(2,4) = -(0);  // 0
    B(2,5) = -(1);  // -1
    B(2,7) = -(1);  // -1
    B(2,8) = -(0);  // 0

    // Now compute local stiffness: K_local = AREA * B^T * D * B
    Eigen::Matrix<double, 9, 9> K_local = AREA * B.transpose() * D * B;

    // Transform to global coordinates: K_global = L^T * K_local * L
    Eigen::Matrix<double, 9, 9> K_global = L.transpose() * K_local * L;

    return K_global;
}

int main() {
    // Compute the stiffness matrix
    Eigen::Matrix<double, 9, 9> K = compute_ktrplt_stiffness();

    // Print as JSON
    print_json_matrix(K);

    return 0;
}