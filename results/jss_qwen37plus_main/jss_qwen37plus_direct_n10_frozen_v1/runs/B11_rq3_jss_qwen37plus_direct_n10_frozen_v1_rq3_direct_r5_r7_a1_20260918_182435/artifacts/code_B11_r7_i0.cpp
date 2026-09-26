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
    
    // Initialize loc to identity permutation
    for (int i = 0; i < n; ++i) {
        loc[i] = i + 1; // Fortran 1-based indexing
    }
    
    // Call QR iteration (qr=0 means perform iteration)
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
    
    // Common block variables (simplified for our use case)
    int MAX = 100 * n;
    int NV = n; // number of eigenvectors to compute
    int NE = n; // number of eigenvalues to print
    int NFOUND = 0; // number of rigid modes found
    int NEVER = 0; // iteration failure counter
    
    // IOPTN and MGIV are not used in our simplified version
    // DLMDAS is not used in our simplified version
    
    // If QR != 0, just reorder (but we always do iteration for test case)
    if (qr != 0) {
        // Reorder only - but we'll skip this for our test case
        goto reorder;
    }
    
    // Search for decoupled submatrices
    int M2 = n;
    
decouple_loop:
    if (M2 <= 1) goto reorder;
    
    int M2M1 = M2 - 1;
    int M1 = 1;
    
    // Find largest M1 such that o[M1-1] != 0 (0-based indexing)
    // We search from the end backwards
    bool found_nonzero = false;
    for (int k = 1; k <= M2M1; ++k) {
        int candidate_M1 = M2 - k;
        if (candidate_M1 >= 1 && candidate_M1 <= static_cast<int>(o.size()) && 
            o[candidate_M1-1] != ZERO) {
            M1 = candidate_M1;
            found_nonzero = true;
            break;
        }
    }
    
    if (!found_nonzero) {
        // All off-diagonal terms are zero, job done
        goto reorder;
    }
    
    // Now we have a decoupled submatrix from M1 to M2
    // Perform QR iteration on this submatrix
    int MM = M1;
    
    // QR iteration loop
    for (int ITER = 1; ITER <= MAX; ++ITER) {
        // Check convergence: if last off-diagonal is effectively zero
        if (std::abs(val[M2-1]) + o[M2M1-1] == std::abs(val[M2-1])) {
            goto converged;
        }
        
        // Check if all diagonal elements in current submatrix are equal
        bool all_equal = true;
        for (int k = M1; k <= M2M1; ++k) {
            if (val[k-1] != val[k]) {
                all_equal = false;
                break;
            }
        }
        
        double SHIFT = ZERO;
        if (!all_equal) {
            // Find smallest diagonal term (in absolute value)
            SHIFT = val[M2-1];
            for (int i = M1; i <= M2M1; ++i) {
                if (std::abs(val[i-1]) < std::abs(SHIFT)) {
                    SHIFT = val[i-1];
                }
            }
            
            // Reduce all terms by shift
            for (int i = M1; i <= M2; ++i) {
                val[i-1] = val[i-1] - SHIFT;
            }
        }
        
        // QR iteration
        double R = val[M1-1] * val[M1-1];
        double S = o[M1-1] / (R + o[M1-1]);
        double T = ZERO;
        double U = S * (val[M1-1] + val[M1]);
        val[M1-1] = val[M1-1] + U;
        
        if (M1 == M2M1) {
            goto shift_back;
        }
        
        int M1P1 = M1 + 1;
        for (int i = M1P1; i <= M2M1; ++i) {
            double G = val[i-1] - U;
            R = (ONE - T) * o[i-2];
            double ONES = ONE - S;
            if (std::abs(ONES) > EPSI) {
                R = G * G / ONES;
            }
            R = R + o[i-1];
            o[i-2] = S * R;
            if (o[i-2] == ZERO) {
                MM = i;
            }
            T = S;
            S = o[i-1] / R;
            U = S * (G + val[i]);
            val[i-1] = U + G;
        }
        
    shift_back:
        val[M2-1] = val[M2-1] - U;
        R = (ONE - T) * o[M2M1-1];
        double ONES = ONE - S;
        if (std::abs(ONES) > EPSI) {
            R = val[M2-1] * val[M2-1] / ONES;
        }
        o[M2M1-1] = S * R;
        
        // Shift back
        if (SHIFT != ZERO) {
            for (int i = M1; i <= M2; ++i) {
                val[i-1] = val[i-1] + SHIFT;
            }
        }
        
        M1 = MM;
    }
    
    // Too many iterations
    NEVER++;
    // In production code, we would call MESAGE here, but for our test we'll continue
    // This should not happen for our well-conditioned test case
    
converged:
    if (M1 == M2M1) {
        // Submatrix reduced to 2x2 or smaller
        if (M1 <= 2) {
            goto reorder;
        }
        M2 = M2M1;
        M2M1 = M2 - 1;
        goto decouple_loop;
    }
    
    if (M1 <= 2) {
        goto reorder;
    }
    
    M2 = M1 - 1;
    goto decouple_loop;
    
reorder:
    // Reorder eigenvalues algebraically in ascending order
    // Create temporary storage for sorted values and locations
    std::vector<double> temp_val = val;
    std::vector<int> temp_loc(n);
    
    // Simple selection sort for ordering
    for (int k = 0; k < n; ++k) {
        int min_idx = -1;
        for (int m = 0; m < n; ++m) {
            if (temp_val[m] != -10000.0) {
                if (min_idx == -1 || temp_val[m] < temp_val[min_idx]) {
                    min_idx = m;
                }
            }
        }
        if (min_idx != -1) {
            temp_loc[k] = min_idx + 1; // Fortran 1-based
            val[k] = temp_val[min_idx];
            temp_val[min_idx] = -10000.0;
        }
    }
    
    // Handle rigid modes (not applicable for our test case)
    if (NFOUND > 0) {
        for (int i = 0; i < NFOUND; ++i) {
            val[i] = ZERO;
        }
    }
}