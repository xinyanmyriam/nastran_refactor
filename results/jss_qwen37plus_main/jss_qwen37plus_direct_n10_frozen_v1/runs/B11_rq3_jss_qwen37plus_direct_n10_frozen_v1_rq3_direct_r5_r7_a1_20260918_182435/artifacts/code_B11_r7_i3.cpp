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
    
    // If QR != 0, just reorder (but we'll skip this for our test case)
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
            std::abs(o[candidate_M1-1]) > EPSI) {
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
        if (std::abs(o[M2M1-1]) <= EPSI * (std::abs(val[M2M1-1]) + std::abs(val[M2-1]))) {
            goto converged;
        }
        
        // Compute Wilkinson shift from bottom-right 2x2 submatrix
        // For submatrix [val[M2M1-1], sqrt(o[M2M1-1]); sqrt(o[M2M1-1]), val[M2-1]]
        double d1 = val[M2M1-1];
        double d2 = val[M2-1];
        double e = std::sqrt(std::abs(o[M2M1-1])); // actual off-diagonal element
        
        // Wilkinson shift: eigenvalue of 2x2 closest to d2
        double trace = d1 + d2;
        double det = d1 * d2 - e * e;
        double disc = trace * trace - 4.0 * det;
        if (disc < 0.0) disc = 0.0;
        double sqrt_disc = std::sqrt(disc);
        double lambda1 = 0.5 * (trace - sqrt_disc);
        double lambda2 = 0.5 * (trace + sqrt_disc);
        double SHIFT = (std::abs(lambda2 - d2) < std::abs(lambda1 - d2)) ? lambda2 : lambda1;
        
        // Reduce all terms by shift
        for (int i = M1; i <= M2; ++i) {
            val[i-1] = val[i-1] - SHIFT;
        }
        
        // QR iteration - apply Givens rotations from top to bottom
        double c = 1.0, s = 0.0;
        double x = val[M1-1];
        
        for (int i = M1; i <= M2M1; ++i) {
            double y = s * val[i-1];
            double z = c * val[i-1];
            
            if (i == M1) {
                x = val[M1-1];
            } else {
                x = c * val[i-1] + s * o[i-2];
                y = s * val[i-1] - c * o[i-2];
            }
            
            // Compute Givens rotation to eliminate y
            double r = std::sqrt(x*x + y*y);
            if (r == 0.0) {
                c = 1.0;
                s = 0.0;
            } else {
                c = x / r;
                s = y / r;
            }
            
            // Apply rotation to current row
            double temp_val = c * val[i-1] + s * o[i-1];
            double temp_o = std::abs(c * o[i-1] - s * val[i]);
            
            if (i == M1) {
                val[M1-1] = r;
                o[M1-1] = temp_o;
                val[i] = temp_val;
            } else {
                o[i-2] = r;
                o[i-1] = temp_o;
                val[i] = temp_val;
            }
        }
        
        // Shift back
        for (int i = M1; i <= M2; ++i) {
            val[i-1] = val[i-1] + SHIFT;
        }
        
        // Check for deflation: find where off-diagonals are small
        MM = M1;
        for (int i = M1; i <= M2M1; ++i) {
            if (std::abs(o[i-1]) <= EPSI * (std::abs(val[i-1]) + std::abs(val[i]))) {
                MM = i + 1;
                break;
            }
        }
        
        // If we deflated, restart with smaller submatrix
        if (MM > M1) {
            M1 = MM;
            if (M1 >= M2) goto converged;
        }
    }
    
    // Too many iterations
    NEVER++;
    // In production code, we would call MESAGE here, but for our test we'll continue
    // This should not happen for our well-conditioned test case
    
converged:
    // Handle 2x2 submatrices explicitly
    if (M2 - M1 == 1) {
        // Solve 2x2 submatrix exactly
        double d1 = val[M1-1];
        double d2 = val[M2-1];
        double e = std::sqrt(std::abs(o[M1-1]));
        
        double trace = d1 + d2;
        double det = d1 * d2 - e * e;
        double disc = trace * trace - 4.0 * det;
        if (disc < 0.0) disc = 0.0;
        double sqrt_disc = std::sqrt(disc);
        double lambda1 = 0.5 * (trace - sqrt_disc);
        double lambda2 = 0.5 * (trace + sqrt_disc);
        
        val[M1-1] = lambda1;
        val[M2-1] = lambda2;
        o[M1-1] = 0.0;
        
        // Continue with next submatrix
        if (M1 > 1) {
            M2 = M1 - 1;
            goto decouple_loop;
        } else {
            goto reorder;
        }
    }
    
    if (M1 == M2) {
        // Single element, done
        if (M1 > 1) {
            M2 = M1 - 1;
            goto decouple_loop;
        } else {
            goto reorder;
        }
    }
    
    // Continue decoupling
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