#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Forward declaration for the main QR iteration function
void qriter(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr);

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    const int n = 5;
    std::vector<double> val = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> o = {1.0, 1.0, 1.0, 1.0}; // squares of off-diagonals: (-1)^2 = 1
    std::vector<int> loc(n);
    
    // Initialize loc to 0-based indices
    for (int i = 0; i < n; ++i) {
        loc[i] = i;
    }
    
    // Call QR iteration (qr = 0 means perform iteration)
    qriter(val, o, loc, 0);
    
    // Sort eigenvalues in ascending order
    std::sort(val.begin(), val.end());
    
    // Output JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < n; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << std::fixed << std::setprecision(12) << val[i];
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}

void qriter(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr) {
    const int n = static_cast<int>(val.size());
    if (n == 0) return;
    
    // Constants
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE = 1.0;
    
    // Set up common variables (simplified from Fortran)
    int max_iter = 100 * n;
    int nv = n; // number of eigenvectors to compute
    int ne = n; // number of eigenvalues to print
    int nfound = 0; // number of rigid modes found
    int ioptn = 0; // not MGIV method
    double dlmdas = 0.0; // not used in standard case
    
    // If qr != 0, just reorder (but we'll do full iteration for test case)
    if (qr != 0) {
        // Just sort the eigenvalues
        std::vector<double> temp = val;
        std::vector<int> temp_loc = loc;
        
        // Simple selection sort to match Fortran logic
        for (int k = 0; k < n; ++k) {
            int m = 0;
            for (int i = 0; i < n; ++i) {
                if (temp[i] != -10000.0) {
                    m = i;
                    break;
                }
            }
            if (m == 0 && temp[0] == -10000.0) break;
            
            for (int i = m + 1; i < n; ++i) {
                if (temp[i] != -10000.0 && temp[m] > temp[i]) {
                    m = i;
                }
            }
            
            o[k] = temp[m];
            temp[m] = -10000.0;
            loc[k] = m;
        }
        
        for (int i = 0; i < n; ++i) {
            val[i] = o[i];
        }
        
        if (nfound > 0) {
            for (int i = 0; i < nfound; ++i) {
                val[i] = ZERO;
            }
        }
        return;
    }
    
    // Main QR iteration loop
    int m1 = 0;
    int m2 = n - 1;
    int never = 0;
    
    // Search for decoupled submatrices
    while (m1 <= m2) {
        // Find largest m2 such that o[m2-1] != 0 (0-based indexing)
        int m2m1 = m2;
        bool found_nonzero = false;
        for (int k = 0; k < m2m1; ++k) {
            int idx = m2m1 - 1 - k; // start from m2-1 down to 0
            if (idx >= 0 && idx < static_cast<int>(o.size()) && o[idx] != ZERO) {
                m2m1 = idx;
                found_nonzero = true;
                break;
            }
        }
        
        if (!found_nonzero) {
            // All off-diagonals are zero in this block
            m1 = m2 + 1;
            continue;
        }
        
        // Now find smallest m1 such that o[m1-1] != 0 (for the block)
        int mm = m2m1 + 1;
        for (int k = m2m1; k > 0; --k) {
            if (k - 1 >= 0 && k - 1 < static_cast<int>(o.size()) && o[k-1] != ZERO) {
                mm = k;
                break;
            }
        }
        m1 = mm;
        
        // QR iteration for submatrix from m1 to m2 (0-based)
        int iter_count = 0;
        bool converged = false;
        
        while (iter_count < max_iter && !converged) {
            iter_count++;
            
            // Check convergence: if |val[m2]| + o[m2m1] == |val[m2]| then o[m2m1] is effectively zero
            if (std::abs(val[m2]) + (m2m1 >= 0 && m2m1 < static_cast<int>(o.size()) ? o[m2m1] : ZERO) == std::abs(val[m2])) {
                converged = true;
                continue;
            }
            
            // Check if all diagonal elements in current block are equal
            bool all_equal = true;
            for (int k = m1; k < m2; ++k) {
                if (val[k] != val[k+1]) {
                    all_equal = false;
                    break;
                }
            }
            
            double shift = ZERO;
            if (!all_equal) {
                // Find smallest diagonal term (in absolute value)
                shift = val[m2];
                for (int i = m1; i <= m2; ++i) {
                    if (std::abs(val[i]) < std::abs(shift)) {
                        shift = val[i];
                    }
                }
                
                // Reduce all terms by shift
                for (int i = m1; i <= m2; ++i) {
                    val[i] -= shift;
                }
            }
            
            // QR iteration
            double r = val[m1] * val[m1];
            double s = (m1 < static_cast<int>(o.size()) ? o[m1] : ZERO) / (r + (m1 < static_cast<int>(o.size()) ? o[m1] : ZERO));
            double t = ZERO;
            double u = s * (val[m1] + (m1 + 1 <= m2 ? val[m1 + 1] : ZERO));
            val[m1] += u;
            
            if (m1 == m2m1) {
                // Single element block
                if (shift != ZERO) {
                    val[m1] += shift;
                }
                continue;
            }
            
            // Main QR loop
            int m1p1 = m1 + 1;
            for (int i = m1p1; i <= m2m1; ++i) {
                double g = val[i] - u;
                r = (ONE - t) * (i - 1 >= 0 && i - 1 < static_cast<int>(o.size()) ? o[i-1] : ZERO);
                double ones = ONE - s;
                if (std::abs(ones) > EPSI) {
                    r = g * g / ones;
                }
                r += (i < static_cast<int>(o.size()) ? o[i] : ZERO);
                if (i - 1 >= 0 && i - 1 < static_cast<int>(o.size())) {
                    o[i-1] = s * r;
                    if (o[i-1] == ZERO) {
                        mm = i;
                    }
                }
                t = s;
                s = (i < static_cast<int>(o.size()) ? o[i] : ZERO) / r;
                u = s * (g + (i + 1 <= m2 ? val[i + 1] : ZERO));
                val[i] = u + g;
            }
            
            val[m2] -= u;
            r = (ONE - t) * (m2m1 >= 0 && m2m1 < static_cast<int>(o.size()) ? o[m2m1] : ZERO);
            double ones = ONE - s;  // Declare ones here for the second usage
            if (std::abs(ones) > EPSI) {
                r = val[m2] * val[m2] / ones;
            }
            if (m2m1 >= 0 && m2m1 < static_cast<int>(o.size())) {
                o[m2m1] = s * r;
            }
            
            // Shift back
            if (shift != ZERO) {
                for (int i = m1; i <= m2; ++i) {
                    val[i] += shift;
                }
            }
            
            m1 = mm;
        }
        
        if (!converged) {
            never++;
        }
        
        // Update bounds for next submatrix
        if (m1 == m2m1) {
            m2 = m2m1 - 1;
        } else {
            m2 = m2m1;
            m2m1 = m2 - 1;
        }
        
        if (m1 <= 1) {
            break;
        }
        m2 = m1 - 1;
        m1 = 0;
    }
    
    // Reorder eigenvalues algebraically in ascending order
    // Using selection sort as in Fortran
    std::vector<double> temp_val = val;
    std::vector<int> temp_loc = loc;
    
    for (int k = 0; k < n; ++k) {
        int m = -1;
        // Find first non-marked element
        for (int i = 0; i < n; ++i) {
            if (temp_val[i] != -10000.0) {
                m = i;
                break;
            }
        }
        if (m == -1) break;
        
        // Find minimum in remaining elements
        for (int i = m + 1; i < n; ++i) {
            if (temp_val[i] != -10000.0 && temp_val[m] > temp_val[i]) {
                m = i;
            }
        }
        
        o[k] = temp_val[m];
        temp_val[m] = -10000.0;
        loc[k] = m;
    }
    
    // Copy back
    for (int i = 0; i < n; ++i) {
        val[i] = o[i];
    }
    
    // Handle rigid modes
    if (nfound > 0) {
        for (int i = 0; i < nfound; ++i) {
            val[i] = ZERO;
        }
    }
}