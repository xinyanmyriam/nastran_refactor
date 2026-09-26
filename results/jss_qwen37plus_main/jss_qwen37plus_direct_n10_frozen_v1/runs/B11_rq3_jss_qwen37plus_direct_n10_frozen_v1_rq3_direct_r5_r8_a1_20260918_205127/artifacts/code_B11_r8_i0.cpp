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
// Output: eigenvalues in sorted ascending order
std::vector<double> qriter(const std::vector<double>& d_in, const std::vector<double>& e_in) {
    int n = static_cast<int>(d_in.size());
    if (n == 0) return {};
    if (n == 1) return {d_in[0]};
    
    // Copy input to working arrays
    std::vector<double> d = d_in;
    std::vector<double> e = e_in; // e[i] is the off-diagonal between d[i] and d[i+1]
    
    // Initialize location array (for reordering)
    std::vector<int> loc(n);
    for (int i = 0; i < n; ++i) {
        loc[i] = i + 1; // 1-based indexing as in Fortran
    }
    
    // Parameters from Fortran
    int max_iter = 100 * n;
    int never = 0;
    int m1 = 0, m2 = n - 1; // 0-based indices
    
    // Main iteration loop
    while (m1 <= m2) {
        // Find largest decoupled submatrix starting from m2
        int m2m1 = m2;
        if (m2 > 0) {
            // Search backwards from m2-1 down to m1 for first non-zero e[k]
            bool found_nonzero = false;
            for (int k = m2 - 1; k >= m1; --k) {
                if (std::abs(e[k]) > EPSI) {
                    m2m1 = k;
                    found_nonzero = true;
                    break;
                }
            }
            if (!found_nonzero) {
                // All off-diagonals are zero in [m1, m2], so we're done for this block
                m1 = m2 + 1;
                continue;
            }
        }
        
        // Now we have a submatrix from m1 to m2 (inclusive)
        // with non-zero e[m2m1] (if m2m1 < m2) or the whole block is coupled
        
        // QR iteration for this submatrix
        int iter;
        for (iter = 0; iter < max_iter; ++iter) {
            // Check convergence: if |e[m2m1]| is negligible compared to |d[m2]|, we're done
            if (m2m1 < m2 && std::abs(e[m2m1]) <= EPSI * (std::abs(d[m2]) + std::abs(d[m2m1]))) {
                break;
            }
            
            // Check if all diagonal elements in current block are equal
            bool all_equal = true;
            for (int k = m1; k < m2; ++k) {
                if (std::abs(d[k] - d[k+1]) > EPSI) {
                    all_equal = false;
                    break;
                }
            }
            
            double shift = ZERO;
            if (!all_equal) {
                // Find smallest diagonal element in current block for shift
                shift = d[m2];
                for (int i = m1; i <= m2; ++i) {
                    if (std::abs(d[i]) < std::abs(shift)) {
                        shift = d[i];
                    }
                }
                
                // Subtract shift from diagonal elements
                for (int i = m1; i <= m2; ++i) {
                    d[i] -= shift;
                }
            }
            
            // QR step
            double r = d[m1] * d[m1];
            double s = (m2m1 >= m1) ? e[m1] / (r + e[m1]) : ZERO;
            double t = ZERO;
            double u = s * (d[m1] + d[m1+1]);
            d[m1] += u;
            
            int mm = m1;
            
            if (m1 < m2m1) {
                // Apply Givens rotations through the submatrix
                for (int i = m1 + 1; i <= m2m1; ++i) {
                    double g = d[i] - u;
                    r = (ONE - t) * e[i-1];
                    double ones = ONE - s;
                    if (std::abs(ones) > EPSI) {
                        r = g * g / ones;
                    }
                    r += e[i];
                    e[i-1] = s * r;
                    if (std::abs(e[i-1]) <= EPSI) {
                        mm = i;
                    }
                    t = s;
                    s = (std::abs(r) > EPSI) ? e[i] / r : ZERO;
                    u = s * (g + d[i+1]);
                    d[i] = u + g;
                }
            }
            
            // Last row update
            if (m2m1 >= m1) {
                d[m2] -= u;
                r = (ONE - t) * e[m2m1];
                double ones = ONE - s;
                if (std::abs(ones) > EPSI) {
                    r = d[m2] * d[m2] / ones;
                }
                e[m2m1] = s * r;
            }
            
            // Add shift back
            if (std::abs(shift) > EPSI) {
                for (int i = m1; i <= m2; ++i) {
                    d[i] += shift;
                }
            }
            
            // Update m1 for next iteration
            m1 = mm;
            
            // Check if we've converged for this submatrix
            if (m1 >= m2m1) {
                break;
            }
        }
        
        // If we didn't converge, adjust m2 and continue
        if (iter >= max_iter) {
            never++;
            // In practice, we'd issue a warning, but for this test we continue
        }
        
        // Move to next submatrix
        if (m1 == m2m1) {
            // This submatrix is done, move to next
            if (m1 < m2) {
                m2 = m2m1 - 1;
            } else {
                m1 = m2 + 1;
            }
        } else {
            m2 = m2m1 - 1;
        }
    }
    
    // Reorder eigenvalues in ascending order (Fortran's bubble sort equivalent)
    std::vector<double> sorted_eigenvalues = d;
    std::vector<int> temp_loc(n);
    for (int i = 0; i < n; ++i) {
        temp_loc[i] = i + 1;
    }
    
    // Simple selection sort to match Fortran behavior
    for (int k = 0; k < n; ++k) {
        int min_idx = k;
        for (int m = k; m < n; ++m) {
            if (sorted_eigenvalues[m] < sorted_eigenvalues[min_idx]) {
                min_idx = m;
            }
        }
        // Swap
        std::swap(sorted_eigenvalues[k], sorted_eigenvalues[min_idx]);
        std::swap(temp_loc[k], temp_loc[min_idx]);
    }
    
    return sorted_eigenvalues;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    std::vector<double> diagonal = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> off_diagonal = {-1.0, -1.0, -1.0, -1.0}; // size n-1
    
    // Compute eigenvalues using QR iteration
    std::vector<double> eigenvalues = qriter(diagonal, off_diagonal);
    
    // Sort the eigenvalues (though our algorithm should return them sorted)
    std::sort(eigenvalues.begin(), eigenvalues.end());
    
    // Generate JSON output
    std::ostringstream json;
    json << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < static_cast<int>(eigenvalues.size()); ++i) {
        if (i > 0) json << ",";
        json << std::fixed << std::setprecision(10) << eigenvalues[i];
    }
    json << "]}";
    
    std::cout << json.str() << std::endl;
    
    return 0;
}