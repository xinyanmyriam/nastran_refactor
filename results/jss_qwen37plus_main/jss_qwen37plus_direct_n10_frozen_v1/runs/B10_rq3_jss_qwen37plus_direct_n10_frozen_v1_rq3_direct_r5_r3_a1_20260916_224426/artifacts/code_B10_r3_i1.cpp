#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <cctype>

// Safe JSON number formatting: no scientific notation, minimal trailing zeros
std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(15) << x;
    std::string s = oss.str();
    
    // Remove trailing zeros
    size_t dot_pos = s.find('.');
    if (dot_pos != std::string::npos) {
        size_t last_non_zero = s.find_last_not_of('0');
        if (last_non_zero != std::string::npos && last_non_zero > dot_pos) {
            s.erase(last_non_zero + 1);
        } else if (last_non_zero == dot_pos) {
            s.erase(dot_pos); // remove "."
        }
    }
    return s;
}

std::string to_json_string(const std::vector<double>& vec) {
    if (vec.empty()) return "[]";
    std::string result = "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) result += ",";
        result += format_double(vec[i]);
    }
    result += "]";
    return result;
}

std::string to_json_matrix(const std::vector<std::vector<double>>& mat) {
    if (mat.empty()) return "[]";
    std::string result = "[";
    for (size_t i = 0; i < mat.size(); ++i) {
        if (i > 0) result += ",";
        result += to_json_string(mat[i]);
    }
    result += "]";
    return result;
}

// INVERD implementation in C++17
// Full pivoting Gaussian elimination for determinant and inverse/solution
void INVERD(int ndim, Eigen::MatrixXd& A, int n, Eigen::MatrixXd& B, int m, 
            double& determinant, int& ising, std::vector<std::vector<int>>& index) {
    // Initialize
    const double EPSI = 1.0e-36;
    
    determinant = 1.0;
    ising = 1; // assume nonsingular until proven otherwise
    
    // Initialize index(n, 3) to zeros: index[i][0]=row, [1]=col, [2]=used flag
    index.clear();
    index.resize(n, std::vector<int>(3, 0));
    
    // Main loop over pivot steps (0 to n-1)
    for (int i = 0; i < n; ++i) {
        // Search for pivot: max |A(j,k)| where row j and col k are unused
        double amax = 0.0;
        int irow = -1;
        int icolum = -1;
        
        for (int j = 0; j < n; ++j) {
            if (index[j][2] == 1) continue; // row j already used
            for (int k = 0; k < n; ++k) {
                if (index[k][2] == 1) continue; // col k already used
                double abs_val = std::abs(A(j, k));
                if (abs_val > amax) {
                    amax = abs_val;
                    irow = j;
                    icolum = k;
                }
            }
        }
        
        // If no pivot found, matrix is singular
        if (irow == -1 || icolum == -1) {
            ising = 2;
            return;
        }
        
        // Mark this column as used (we use column index for tracking)
        index[icolum][2] = 1;
        // Store pivot row and column for later column reordering
        index[i][0] = irow;
        index[i][1] = icolum;
        
        // Interchange rows irow and icolum to bring pivot to (icolum, icolum)
        if (irow != icolum) {
            determinant = -determinant;
            A.row(irow).swap(A.row(icolum));
            if (m > 0) {
                B.row(irow).swap(B.row(icolum));
            }
        }
        
        // Pivot element
        double pivot = A(icolum, icolum);
        determinant *= pivot;
        
        if (std::abs(pivot) < EPSI) {
            ising = 2;
            return;
        }
        
        // Normalize pivot row
        A.row(icolum) /= pivot;
        if (m > 0) {
            B.row(icolum) /= pivot;
        }
        
        // Eliminate other rows
        for (int l1 = 0; l1 < n; ++l1) {
            if (l1 == icolum) continue;
            double t = A(l1, icolum);
            if (std::abs(t) < EPSI) continue;
            
            A.row(l1) -= t * A.row(icolum);
            if (m > 0) {
                B.row(l1) -= t * B.row(icolum);
            }
        }
    }
    
    // Apply column interchanges (reverse order) to reconstruct inverse
    // We stored (irow, icolum) for step i; now undo column swaps from last to first
    for (int i = n-1; i >= 0; --i) {
        int jrow = index[i][0];
        int jcolum = index[i][1];
        if (jrow != jcolum) {
            A.col(jrow).swap(A.col(jcolum));
        }
    }
    
    // Final singularity check: ensure all columns were used
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
    
    // A now contains the inverse of K
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
    std::cout << ",\"determinant\":" << format_double(determinant);
    std::cout << ",\"inverse\":" << to_json_matrix(inverse_matrix) << "}";
    std::cout << std::endl;
    
    return 0;
}