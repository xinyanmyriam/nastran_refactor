#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

// Define basic types for compatibility
using Real = double;

// Simple 3D vector class
class Vector3d {
public:
    Real x, y, z;
    
    Vector3d(Real x_ = 0.0, Real y_ = 0.0, Real z_ = 0.0) : x(x_), y(y_), z(z_) {}
    
    // Vector operations needed for the stiffness computation
    Vector3d operator-(const Vector3d& other) const {
        return Vector3d(x - other.x, y - other.y, z - other.z);
    }
    
    Real norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};

// Simple 12x12 matrix class
class Matrix12d {
private:
    Real data[12][12];
    
public:
    Matrix12d() {
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                data[i][j] = 0.0;
            }
        }
    }
    
    Real& operator()(int i, int j) {
        return data[i][j];
    }
    
    const Real& operator()(int i, int j) const {
        return data[i][j];
    }
};

// Mock MAT function to provide material properties
void MAT(int id, Real& E, Real& G, Real& nu, Real& rho, Real& alpha) {
    // Typical steel properties
    E = 200.0e9;    // Pa
    G = 77.0e9;     // Pa
    nu = 0.3;
    rho = 7850.0;   // kg/m^3
    alpha = 12.0e-6; // 1/K
}

// Helper function to format double for JSON output
std::string format_double(Real value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << value;
    std::string result = oss.str();
    
    // Remove trailing zeros and potential decimal point
    result.erase(result.find_last_not_of('0') + 1, std::string::npos);
    if (result.back() == '.') {
        result.pop_back();
    }
    
    return result;
}

// Simplified CROD stiffness matrix computation
// For a simple rod element: axial stiffness only (no torsion in basic CROD)
Matrix12d compute_crod_stiffness(const Vector3d& nodeA, const Vector3d& nodeB, 
                                Real E, Real A, Real G, Real J) {
    Matrix12d K;
    
    // Compute element length
    Vector3d diff = nodeB - nodeA;
    Real L = diff.norm();
    
    if (L <= 0.0) {
        throw std::invalid_argument("Element length must be positive");
    }
    
    // Axial stiffness
    Real k_axial = E * A / L;
    
    // Torsional stiffness
    Real k_torsion = G * J / L;
    
    // For CROD element with 2 nodes (each with 6 DOF: tx,ty,tz,rx,ry,rz):
    // Node A: DOFs 0,1,2,3,4,5
    // Node B: DOFs 6,7,8,9,10,11
    //
    // Axial stiffness (x-direction) affects DOF 0 (tx_A) and DOF 6 (tx_B)
    K(0,0) = k_axial;
    K(0,6) = -k_axial;
    K(6,0) = -k_axial;
    K(6,6) = k_axial;
    
    // Torsional stiffness (x-direction rotation) affects DOF 3 (rx_A) and DOF 9 (rx_B)
    K(3,3) = k_torsion;
    K(3,9) = -k_torsion;
    K(9,3) = -k_torsion;
    K(9,9) = k_torsion;
    
    return K;
}

int main() {
    try {
        // Test case data
        Vector3d nodeA(0.0, 0.0, 0.0);
        Vector3d nodeB(2.0, 0.0, 0.0);
        Real E, G, nu, rho, alpha;
        MAT(1, E, G, nu, rho, alpha); // Get E and G from MAT
        Real A = 0.01;         // m^2
        Real J = 5.0e-6;       // m^4
        
        // Compute stiffness matrix
        Matrix12d K = compute_crod_stiffness(nodeA, nodeB, E, A, G, J);
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 12; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << "[";
            for (int j = 0; j < 12; ++j) {
                if (j > 0) std::cout << ",";
                std::cout << format_double(K(i,j));
            }
            std::cout << "]";
        }
        std::cout << "]}" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}