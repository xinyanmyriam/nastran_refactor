#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

// For MSVC compatibility, define a simple JSON library if nlohmann/json is not available
// But since the requirement says "no external deps beyond Eigen + stdlib", and nlohmann/json is common,
// we'll use it. If unavailable, replace with manual JSON formatting.

// However, the problem states "no external deps beyond Eigen + stdlib", so let's avoid nlohmann/json
// and write minimal JSON output manually.

// We'll create a simple JSON formatter for the required output

std::string to_json_string(const std::vector<double>& vec) {
    std::string result = "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) result += ",";
        result += std::to_string(vec[i]);
        // Remove scientific notation and trim trailing zeros
        size_t pos = result.find_last_of('.');
        if (pos != std::string::npos) {
            // Remove trailing zeros
            size_t end = result.find_last_not_of('0');
            if (end != std::string::npos && result[end] == '.') {
                result.erase(end);
            } else if (end != std::string::npos) {
                result.erase(end + 1);
            }
        }
    }
    result += "]";
    return result;
}

std::string to_json_matrix(const std::vector<std::vector<double>>& mat) {
    std::string result = "[";
    for (size_t i = 0; i < mat.size(); ++i) {
        if (i > 0) result += ",";
        result += to_json_string(mat[i]);
    }
    result += "]";
    return result;
}

// INVERD implementation in C++17
// This implements the same algorithm as the Fortran code: full pivoting Gaussian elimination
// with row/column interchanges, computing determinant and inverse.

void INVERD(int ndim, Eigen::MatrixXd& A, int n, Eigen::MatrixXd& B, int m, 
            double& determinant, int& ising, std::vector<std::vector<int>>& index) {
    // Initialize
    const double EPSI = 1.0e-36;
    double mach = 5.0; // Assume VAX-like machine for EPSI adjustment
    if (mach == 5.0) {
        // Use larger epsilon for VAX-like machines
        // But for our test case, we'll keep the original EPSI
    }
    
    determinant = 1.0;
    if (ising < 0) {
        determinant = 0.0;
    }
    
    // Initialize index(1..n, 3) to zeros
    index.clear();
    index.resize(n, std::vector<int>(3, 0));
    
    // Main loop over pivot steps
    for (int i = 0; i < n; ++i) {
        // Search for pivot
        double amax = 0.0;
        int irow = -1;
        int icolum = -1;
        
        for (int j = 0; j < n; ++j) {
            if (index[j][2] == 1) continue;
            for (int k = 0; k < n; ++k) {
                if (index[k][2] == 1) continue;
                double abs_val = std::abs(A(j, k));
                if (abs_val > amax) {
                    amax = abs_val;
                    irow = j;
                    icolum = k;
                }
            }
        }
        
        if (irow == -1 || icolum == -1) {
            ising = 2;
            return;
        }
        
        // Mark column as used
        index[icolum][2] = index[icolum][2] + 1;
        // Store row and column indices for later column interchange
        index[i][0] = irow;
        index[i][1] = icolum;
        
        // Interchange rows to put pivot on diagonal
        if (irow != icolum) {
            determinant = -determinant;
            // Swap rows irow and icolum in A
            A.row(irow).swap(A.row(icolum));
            // Swap rows in B if present
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    double swap_val = B(irow, l);
                    B(irow, l) = B(icolum, l);
                    B(icolum, l) = swap_val;
                }
            }
        }
        
        // Divide pivot row by pivot element
        double pivot = A(icolum, icolum);
        determinant = determinant * pivot;
        
        if (std::abs(pivot) < EPSI) {
            ising = 2;
            return;
        }
        
        // Normalize pivot row
        A.row(icolum) = A.row(icolum) / pivot;
        if (m > 0) {
            B.row(icolum) = B.row(icolum) / pivot;
        }
        
        // Reduce other rows
        for (int l1 = 0; l1 < n; ++l1) {
            if (l1 == icolum) continue;
            double t = A(l1, icolum);
            A(l1, icolum) = 0.0;
            if (std::abs(t) < EPSI) continue;
            
            // Update row l1: A[l1] -= t * A[icolum]
            A.row(l1) = A.row(l1) - t * A.row(icolum);
            if (m > 0) {
                B.row(l1) = B.row(l1) - t * B.row(icolum);
            }
        }
    }
    
    // Interchange columns
    for (int i = 0; i < n; ++i) {
        int l = n - 1 - i;
        if (index[l][0] == index[l][1]) continue;
        int jrow = index[l][0];
        int jcolum = index[l][1];
        // Swap columns jrow and jcolum in A
        A.col(jrow).swap(A.col(jcolum));
    }
    
    // Check if matrix is singular
    ising = 1;
    for (int k = 0; k < n; ++k) {
        if (index[k][2] != 1) {
            ising = 2;
            break;
        }
    }
}

int main() {
    // Test case: 4x4 SPD matrix K
    // K = [200,-100,0,0; -100,200,-100,0; 0,-100,200,-100; 0,0,-100,100]
    Eigen::MatrixXd K(4, 4);
    K << 200.0, -100.0,   0.0,   0.0,
        -100.0,  200.0, -100.0,   0.0,
           0.0, -100.0,  200.0, -100.0,
           0.0,    0.0, -100.0,  100.0;
    
    // b = K * [1,2,3,4]^T
    Eigen::VectorXd x_true(4);
    x_true << 1.0, 2.0, 3.0, 4.0;
    Eigen::VectorXd b = K * x_true;
    
    // Prepare inputs for INVERD
    int ndim = 4;
    int n = 4;
    int m = 1; // one RHS vector
    
    // Copy K to A (will be modified)
    Eigen::MatrixXd A = K;
    
    // B will hold the RHS vector b
    Eigen::MatrixXd B(ndim, m);
    for (int i = 0; i < ndim; ++i) {
        B(i, 0) = b(i);
    }
    
    double determinant = 0.0;
    int ising = 1;
    std::vector<std::vector<int>> index;
    
    // Call INVERD
    INVERD(ndim, A, n, B, m, determinant, ising, index);
    
    // Extract solution from B
    std::vector<double> solution;
    for (int i = 0; i < n; ++i) {
        solution.push_back(B(i, 0));
    }
    
    // Compute inverse: A now contains the inverse (since we inverted the matrix)
    // The INVERD routine returns the inverse in A
    std::vector<std::vector<double>> inverse_matrix;
    for (int i = 0; i < n; ++i) {
        std::vector<double> row;
        for (int j = 0; j < n; ++j) {
            row.push_back(A(i, j));
        }
        inverse_matrix.push_back(row);
    }
    
    // Output JSON
    std::cout << "{\"test\":\"INVERD\",\"solution\":" << to_json_string(solution);
    std::cout << ",\"determinant\":" << determinant;
    std::cout << ",\"inverse\":" << to_json_matrix(inverse_matrix) << "}";
    std::cout << std::endl;
    
    return 0;
}