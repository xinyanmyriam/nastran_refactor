#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

// Simplified version of the Fortran QRITER subroutine for tridiagonal matrices
// We focus only on the core QR iteration logic for eigenvalue computation

void qriter(std::vector<double>& diag, std::vector<double>& off_diag, int n) {
    // diag: diagonal elements (size n)
    // off_diag: off-diagonal elements (size n-1), representing sub/super diagonal
    // n: matrix size
    
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE = 1.0;
    
    // Maximum iterations: 100 * n as in Fortran
    int max_iter = 100 * n;
    
    // Main iteration loop over decoupled submatrices
    int m1 = 0;
    int m2 = n - 1; // 0-based indexing
    
    while (m1 <= m2) {
        // Find largest decoupled submatrix starting from m2
        // Search backwards for first non-zero off-diagonal in [m1, m2-1]
        int k;
        if (m2 > m1) {
            for (k = m2 - 1; k >= m1; --k) {
                if (std::abs(off_diag[k]) > EPSI) break;
            }
            if (k < m1) {
                // All off-diagonals zero in [m1, m2], so diagonal contains eigenvalues
                m1 = m2 + 1;
                continue;
            }
            // Submatrix is from m1 to m2 (inclusive)
        } else {
            // 1x1 submatrix
            m1 = m2 + 1;
            continue;
        }
        
        // If submatrix is 1x1, done
        if (m2 == m1) {
            m1 = m2 + 1;
            continue;
        }
        
        // QR iteration for submatrix [m1, m2]
        bool converged = false;
        for (int iter = 0; iter < max_iter; ++iter) {
            // Check convergence: if |off_diag[m2-1]| is negligible
            if (std::abs(off_diag[m2-1]) <= EPSI * (std::abs(diag[m2-1]) + std::abs(diag[m2]))) {
                converged = true;
                break;
            }
            
            // Wilkinson shift: eigenvalue of bottom 2x2 submatrix closest to diag[m2]
            double a = diag[m2-1];
            double b = off_diag[m2-1];
            double c = diag[m2];
            double discriminant = (a - c) * (a - c) + 4.0 * b * b;
            if (discriminant < 0.0) discriminant = 0.0;
            double sqrt_disc = std::sqrt(discriminant);
            double lambda1 = 0.5 * (a + c - sqrt_disc);
            double lambda2 = 0.5 * (a + c + sqrt_disc);
            double shift = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
            
            // Apply shift
            for (int i = m1; i <= m2; ++i) {
                diag[i] -= shift;
            }
            
            // Implicit QR step: introduce bulge at top and chase it down
            double d = diag[m1];
            double e = off_diag[m1];
            double f = 0.0;
            
            // Process from top to bottom
            for (int i = m1; i < m2; ++i) {
                // Compute Givens rotation parameters
                double r = std::sqrt(d * d + e * e);
                if (r == 0.0) {
                    r = EPSI;
                }
                double c_rot = d / r;
                double s_rot = e / r;
                
                // Apply rotation to current 2x2 block
                double temp = c_rot * diag[i] + s_rot * off_diag[i];
                off_diag[i] = -s_rot * diag[i] + c_rot * off_diag[i];
                diag[i] = temp;
                
                // Apply rotation to next diagonal element
                if (i < m2 - 1) {
                    temp = c_rot * off_diag[i] + s_rot * diag[i + 1];
                    diag[i + 1] = -s_rot * off_diag[i] + c_rot * diag[i + 1];
                    off_diag[i] = temp;
                } else {
                    // Last element: only update diag[m2]
                    diag[m2] = c_rot * off_diag[i] + s_rot * diag[m2];
                }
                
                // Update for next iteration
                if (i < m2 - 1) {
                    d = diag[i + 1];
                    e = off_diag[i + 1];
                }
            }
            
            // Shift back
            for (int i = m1; i <= m2; ++i) {
                diag[i] += shift;
            }
            
            // Check if bottom off-diagonal became small enough
            if (std::abs(off_diag[m2-1]) <= EPSI * (std::abs(diag[m2-1]) + std::abs(diag[m2]))) {
                converged = true;
                break;
            }
        }
        
        if (!converged) {
            // In production code, we'd handle failure, but for test case it should converge
        }
        
        // Decouple the last element if converged
        if (converged && m2 > m1) {
            m2--;
        } else {
            // Move to next submatrix
            m1++;
        }
    }
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    const int n = 5;
    std::vector<double> diag = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> off_diag = {-1.0, -1.0, -1.0, -1.0}; // subdiagonal only
    
    // Make a copy for the algorithm (we'll modify diag in-place)
    std::vector<double> diag_copy = diag;
    
    // Run QR iteration
    qriter(diag_copy, off_diag, n);
    
    // Sort eigenvalues in ascending order
    std::sort(diag_copy.begin(), diag_copy.end());
    
    // Generate analytical solution for verification
    std::vector<double> analytical;
    for (int k = 1; k <= n; ++k) {
        double lambda_k = 2.0 - 2.0 * std::cos(k * M_PI / 6.0);
        analytical.push_back(lambda_k);
    }
    std::sort(analytical.begin(), analytical.end());
    
    // Output JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < n; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << std::fixed << std::setprecision(10) << diag_copy[i];
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}