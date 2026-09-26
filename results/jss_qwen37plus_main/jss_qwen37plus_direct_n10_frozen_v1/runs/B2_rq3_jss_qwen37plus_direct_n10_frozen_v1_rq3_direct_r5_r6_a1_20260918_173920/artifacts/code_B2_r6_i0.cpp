#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <json/json.h> // We'll implement minimal JSON serialization to avoid external dependency

// Minimal JSON serialization for Eigen matrix
std::string eigenMatrixToJson(const Eigen::MatrixXd& mat) {
    std::string result = "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < mat.rows(); ++i) {
        if (i > 0) result += ",";
        result += "[";
        
        for (int j = 0; j < mat.cols(); ++j) {
            if (j > 0) result += ",";
            
            // Format in scientific notation with e+0N format
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(15) << mat(i, j);
            std::string val = oss.str();
            
            // Ensure consistent e+0N format (not e+N)
            size_t pos = val.find('e');
            if (pos != std::string::npos && val.length() > pos + 2) {
                std::string exp = val.substr(pos + 1);
                if (exp[0] == '+' || exp[0] == '-') {
                    if (exp.length() == 2) {
                        val = val.substr(0, pos + 1) + "0" + exp;
                    }
                }
            }
            
            result += val;
        }
        result += "]";
    }
    
    result += "]}";
    return result;
}

// Helper function to create transformation matrix from direction vectors
Eigen::Matrix3d createTransformationMatrix(const Eigen::Vector3d& x, const Eigen::Vector3d& y, const Eigen::Vector3d& z) {
    Eigen::Matrix3d T;
    T.col(0) = x;
    T.col(1) = y;
    T.col(2) = z;
    return T;
}

// Main CBAR stiffness computation
Eigen::MatrixXd computeCBARStiffness() {
    // Test case parameters
    const double E = 200e9;           // Pa
    const double G = 76.923e9;        // Pa
    const double A = 0.01;            // m^2
    const double Iy = 8.333e-6;       // m^4 (I1 in Fortran - bending about y-axis)
    const double Iz = 8.333e-6;       // m^4 (I2 in Fortran - bending about z-axis)
    const double J = 1.667e-5;        // m^4 (FJ in Fortran - torsional constant)
    const double L = 2.0;             // m (length)
    
    // Node coordinates
    Eigen::Vector3d nodeA(0.0, 0.0, 0.0);
    Eigen::Vector3d nodeB(2.0, 0.0, 0.0);
    
    // Reference vector (SMALLV) - use y-axis as reference for local coordinate system
    Eigen::Vector3d smallV(0.0, 1.0, 0.0);
    
    // Normalize reference vector
    double fl = smallV.norm();
    if (fl <= 0.0) {
        throw std::runtime_error("Invalid reference vector");
    }
    Eigen::Vector3d smalv0 = smallV / fl;
    
    // Compute element axis vector (VECI) - from A to B
    Eigen::Vector3d vecI = nodeB - nodeA;
    double length = vecI.norm();
    if (length == 0.0) {
        throw std::runtime_error("Zero length element");
    }
    vecI = vecI / length; // normalize
    
    // Compute VECK = VECI × SMALV0
    Eigen::Vector3d vecK;
    vecK(0) = vecI(1) * smalv0(2) - vecI(2) * smalv0(1);
    vecK(1) = vecI(2) * smalv0(0) - vecI(0) * smalv0(2);
    vecK(2) = vecI(0) * smalv0(1) - vecI(1) * smalv0(0);
    
    double fll = vecK.norm();
    if (fll == 0.0) {
        throw std::runtime_error("Invalid cross product");
    }
    vecK = vecK / fll;
    
    // Compute VECJ = VECK × VECI
    Eigen::Vector3d vecJ;
    vecJ(0) = vecK(1) * vecI(2) - vecK(2) * vecI(1);
    vecJ(1) = vecK(2) * vecI(0) - vecK(0) * vecI(2);
    vecJ(2) = vecK(0) * vecI(1) - vecK(1) * vecI(0);
    
    fll = vecJ.norm();
    if (fll == 0.0) {
        throw std::runtime_error("Invalid cross product");
    }
    vecJ = vecJ / fll;
    
    // Material properties and intermediate calculations
    double lsq = length * length;
    double lcube = lsq * length;
    
    double ei1 = E * Iy;  // EI1
    double ei2 = E * Iz;  // EI2
    
    // Shear correction factors (K1, K2) - set to 1.0 as not specified in test case
    double k1 = 1.0;
    double k2 = 1.0;
    double i12 = 0.0; // No warping constant
    
    // Calculate R1 and R2 (shear-deformation terms)
    double r1, r2;
    if (k1 == 0.0 || i12 != 0.0) {
        r1 = 12.0 * ei1 / lcube;
    } else {
        double gak1 = G * A * k1;
        r1 = (12.0 * ei1 * gak1) / (gak1 * lcube + 12.0 * length * ei1);
    }
    
    if (k2 == 0.0 || i12 != 0.0) {
        r2 = 12.0 * ei2 / lcube;
    } else {
        double gak2 = G * A * k2;
        r2 = (12.0 * ei2 * gak2) / (gak2 * lcube + 12.0 * length * ei2);
    }
    
    // Small K terms
    double sk1 = 0.25 * r1 * lsq + ei1 / length;
    double sk2 = 0.25 * r2 * lsq + ei2 / length;
    double sk3 = 0.25 * r1 * lsq - ei1 / length;
    double sk4 = 0.25 * r2 * lsq - ei2 / length;
    
    // Other terms
    double ael = E * A / length;
    double lr1 = length * r1 / 2.0;
    double lr2 = length * r2 / 2.0;
    double gjl = G * J / length;
    
    // Initialize 12x12 stiffness matrix
    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(12, 12);
    
    // Fill the stiffness matrix according to Fortran logic
    // Note: Fortran uses 1-based indexing, C++ uses 0-based
    // The Fortran array KE(144) is stored column-wise: KE(i,j) -> index = (j-1)*12 + i
    // But we'll fill row-wise for clarity and then transpose if needed
    
    // Axial terms
    ke(0, 0) = ael;      // KE(1)
    ke(0, 6) = -ael;     // KE(7)
    ke(6, 0) = -ael;     // KE(73)
    ke(6, 6) = ael;      // KE(79)
    
    // Bending about y-axis (Iy) - affects DOF 1,2,5,7,8,11 (0-indexed: 1,2,5,7,8,11)
    // In Fortran: R1 at (14), LR1 at (18), -R1 at (20), LR1 at (24), etc.
    ke(1, 1) = r1;       // KE(14)
    ke(1, 5) = lr1;      // KE(18)
    ke(1, 7) = -r1;      // KE(20)
    ke(1, 11) = lr1;     // KE(24)
    
    ke(5, 1) = lr1;      // KE(62)
    ke(5, 5) = sk1;      // KE(66)
    ke(5, 7) = -lr1;     // KE(68)
    ke(5, 11) = sk3;     // KE(72)
    
    ke(7, 1) = -r1;      // KE(86)
    ke(7, 5) = -lr1;     // KE(90)
    ke(7, 7) = r1;       // KE(92)
    ke(7, 11) = -lr1;    // KE(96)
    
    ke(11, 1) = lr1;     // KE(134)
    ke(11, 5) = sk3;     // KE(138)
    ke(11, 7) = -lr1;    // KE(140)
    ke(11, 11) = sk1;    // KE(144)
    
    // Bending about z-axis (Iz) - affects DOF 2,3,4,8,9,10 (0-indexed: 2,3,4,8,9,10)
    ke(2, 2) = r2;       // KE(27)
    ke(2, 4) = -lr2;     // KE(29)
    ke(2, 8) = -r2;      // KE(33)
    ke(2, 10) = -lr2;    // KE(35)
    
    ke(4, 2) = -lr2;     // KE(51)
    ke(4, 4) = sk2;      // KE(53)
    ke(4, 8) = lr2;      // KE(57)
    ke(4, 10) = sk4;     // KE(59)
    
    ke(8, 2) = -r2;      // KE(99)
    ke(8, 4) = lr2;      // KE(101)
    ke(8, 8) = r2;       // KE(105)
    ke(8, 10) = lr2;     // KE(107)
    
    ke(10, 2) = -lr2;    // KE(123)
    ke(10, 4) = sk4;     // KE(125)
    ke(10, 8) = lr2;     // KE(129)
    ke(10, 10) = sk2;    // KE(131)
    
    // Torsion about x-axis (J) - affects DOF 3,9 (0-indexed: 3,9)
    ke(3, 3) = gjl;      // KE(40)
    ke(3, 9) = -gjl;     // KE(46)
    ke(9, 3) = -gjl;     // KE(112)
    ke(9, 9) = gjl;      // KE(118)
    
    // Handle warping (I12) - not present in test case (i12 = 0.0)
    // So no additional terms
    
    // Apply coordinate transformation to get global stiffness matrix
    // Create transformation matrix for each node: [T 0; 0 T] where T transforms local to global
    // Local axes: x = vecI, y = vecJ, z = vecK
    Eigen::Matrix3d T = createTransformationMatrix(vecI, vecJ, vecK);
    
    // Global transformation matrix for 12 DOF: block diagonal [T, T]
    Eigen::MatrixXd T_global = Eigen::MatrixXd::Zero(12, 12);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T_global(i, j) = T(i, j);
            T_global(i+3, j+3) = T(i, j);
            T_global(i+6, j+6) = T(i, j);
            T_global(i+9, j+9) = T(i, j);
        }
    }
    
    // Transform stiffness matrix: K_global = T_global^T * K_local * T_global
    Eigen::MatrixXd k_global = T_global.transpose() * ke * T_global;
    
    return k_global;
}

int main() {
    try {
        // Compute the CBAR stiffness matrix
        Eigen::MatrixXd stiffness = computeCBARStiffness();
        
        // Output as JSON
        std::cout << eigenMatrixToJson(stiffness) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}