#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Forward declarations for Eigen types
using Eigen::VectorXd;
using Eigen::MatrixXd;

// Constants
const double EPSI = 1.0e-10;
const double ZERO = 0.0;
const double ONE = 1.0;

// QR iteration for tridiagonal matrix eigenvalues
// Input: diagonal vector 'd' and off-diagonal vector 'e' (size n-1)
// Output: eigenvalues in 'd' (sorted ascending)
void QRITER(VectorXd& d, const VectorXd& e) {
    int n = static_cast<int>(d.size());
    if (n <= 0) return;
    
    // Handle trivial case
    if (n == 1) return;
    
    // Copy off-diagonal squares as needed (Fortran uses O as squares of off-diagonals)
    VectorXd o = e.array().square();
    
    // LOC array for reordering (we'll use it internally)
    std::vector<int> loc(n);
    for (int i = 0; i < n; ++i) {
        loc[i] = i + 1; // 1-based indexing as in Fortran
    }
    
    // Simulate common blocks
    int MAX = 100 * n;
    int NEVER = 0;
    int M1 = 0, M2 = n - 1;
    
    // Main iteration loop - process submatrices
    while (M1 <= M2) {
        // Find decoupled submatrix starting from M2
        int M2M1 = M2;
        bool found_nonzero = false;
        
        // Search backward from M2 to find first nonzero off-diagonal
        for (int k = 1; k <= M2 - M1; ++k) {
            int idx = M2 - k;
            if (idx >= 0 && o(idx) != ZERO) {
                M2M1 = idx;
                found_nonzero = true;
                break;
            }
        }
        
        if (!found_nonzero) {
            // All off-diagonals zero in this block - eigenvalues are just diagonals
            break;
        }
        
        // Set up current submatrix [M1, M2]
        int mm = M1;
        int iter_count = 0;
        
        // QR iteration for current submatrix
        for (iter_count = 0; iter_count < MAX; ++iter_count) {
            // Check convergence: if last off-diagonal is negligible relative to diagonal
            if (std::abs(d(M2)) + o(M2M1) == std::abs(d(M2))) {
                break;
            }
            
            // Check if all diagonals in current submatrix are equal
            bool all_equal = true;
            for (int k = M1; k < M2; ++k) {
                if (d(k) != d(k+1)) {
                    all_equal = false;
                    break;
                }
            }
            
            double shift = ZERO;
            if (!all_equal) {
                // Find smallest diagonal term (in absolute value) for shift
                shift = d(M2);
                for (int i = M1; i < M2; ++i) {
                    if (std::abs(d(i)) < std::abs(shift)) {
                        shift = d(i);
                    }
                }
                
                // Reduce all terms by shift
                for (int i = M1; i <= M2; ++i) {
                    d(i) -= shift;
                }
            }
            
            // QR iteration step
            double r = d(M1) * d(M1);
            double s = o(M1) / (r + o(M1));
            double t = ZERO;
            double u = s * (d(M1) + d(M1+1));
            d(M1) += u;
            
            if (M1 == M2M1) {
                // Special case: 2x2 submatrix
                d(M2) -= u;
                r = (ONE - t) * o(M2M1);
                double ones = ONE - s;
                if (std::abs(ones) > EPSI) {
                    r = d(M2) * d(M2) / ones;
                }
                o(M2M1) = s * r;
            } else {
                // General case
                int M1P1 = M1 + 1;
                for (int i = M1P1; i <= M2M1; ++i) {
                    double g = d(i) - u;
                    r = (ONE - t) * o(i-1);
                    double ones = ONE - s;
                    if (std::abs(ones) > EPSI) {
                        r = g * g / ones;
                    }
                    r = r + o(i);
                    o(i-1) = s * r;
                    if (o(i-1) == ZERO) {
                        mm = i;
                    }
                    t = s;
                    s = o(i) / r;
                    u = s * (g + d(i+1));
                    d(i) = u + g;
                }
                
                d(M2) -= u;
                r = (ONE - t) * o(M2M1);
                double ones = ONE - s;
                if (std::abs(ones) > EPSI) {
                    r = d(M2) * d(M2) / ones;
                }
                o(M2M1) = s * r;
            }
            
            // Shift back if needed
            if (shift != ZERO) {
                for (int i = M1; i <= M2; ++i) {
                    d(i) += shift;
                }
            }
            
            M1 = mm;
        }
        
        if (iter_count >= MAX) {
            // Failed to converge - increment NEVER counter
            NEVER++;
        }
        
        // Update submatrix bounds
        if (M1 == M2M1) {
            // Move to next submatrix
            M2 = M2M1 - 1;
        } else {
            M2 = M2M1;
            M2M1 = M2 - 1;
        }
        
        if (M1 <= 1) {
            break;
        }
        M2 = M1 - 1;
        M1 = 0;
    }
    
    // Reorder eigenvalues algebraically in ascending order
    // Create temporary arrays for sorting
    VectorXd temp_d = d;
    VectorXd temp_o(n);
    
    // Use -10000.0 as sentinel value (as in Fortran)
    const double SENTINEL = -10000.0;
    
    for (int k = 0; k < n; ++k) {
        // Find minimum element not yet selected
        int min_idx = -1;
        for (int m = 0; m < n; ++m) {
            if (temp_d(m) != SENTINEL) {
                if (min_idx == -1 || temp_d(m) < temp_d(min_idx)) {
                    min_idx = m;
                }
            }
        }
        
        if (min_idx != -1) {
            temp_o(k) = temp_d(min_idx);
            temp_d(min_idx) = SENTINEL;
            loc[k] = min_idx + 1; // Convert to 1-based
        }
    }
    
    // Copy sorted eigenvalues back
    d = temp_o;
}

// Analytical solution for test case: lambda_k = 2 - 2*cos(k*pi/6), k=1..5
std::vector<double> analytical_eigenvalues() {
    std::vector<double> result;
    for (int k = 1; k <= 5; ++k) {
        double lambda = 2.0 - 2.0 * std::cos(k * M_PI / 6.0);
        result.push_back(lambda);
    }
    return result;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    int n = 5;
    VectorXd d(n);
    VectorXd e(n-1);
    
    // Initialize diagonal
    for (int i = 0; i < n; ++i) {
        d(i) = 2.0;
    }
    
    // Initialize off-diagonal (subdiagonal, symmetric)
    for (int i = 0; i < n-1; ++i) {
        e(i) = -1.0;
    }
    
    // Call QRITER
    QRITER(d, e);
    
    // Sort eigenvalues in ascending order (QRITER should have done this,
    // but ensure for JSON output)
    std::vector<double> eigenvalues;
    for (int i = 0; i < n; ++i) {
        eigenvalues.push_back(d(i));
    }
    std::sort(eigenvalues.begin(), eigenvalues.end());
    
    // Format JSON output
    std::ostringstream json;
    json << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (i > 0) json << ",";
        json << std::fixed << std::setprecision(12) << eigenvalues[i];
    }
    json << "]}";
    
    std::cout << json.str() << std::endl;
    
    return 0;
}