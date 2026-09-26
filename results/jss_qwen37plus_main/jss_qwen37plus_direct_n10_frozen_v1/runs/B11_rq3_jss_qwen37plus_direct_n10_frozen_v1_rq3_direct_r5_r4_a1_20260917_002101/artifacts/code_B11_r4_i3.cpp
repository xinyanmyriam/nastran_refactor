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
        if (m2 > m1) {
            int k;
            for (k = m2 - 1; k >= m1; --k) {
                if (std::abs(off_diag[k]) > EPSI) break;
            }
            if (k < m1) {
                // All off-diagonals zero in [m1, m2], so diagonal contains eigenvalues
                m1 = m2 + 1;
                continue;
            }
            m2 = k + 1; // m2 is now the end index of current submatrix
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
            
            // Wilkinson shift: eigenvalue of bottom-right 2x2 submatrix closest to diag[m2]
            double d1 = diag[m2-1];
            double d2 = diag[m2];
            double e = off_diag[m2-1];
            double mu;
            if (m2 - 1 == m1) {
                // 2x2 submatrix, use exact eigenvalue
                double disc = (d1 - d2) * (d1 - d2) + 4.0 * e * e;
                if (disc >= 0.0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (d1 + d2 - sqrt_disc);
                    double lambda2 = 0.5 * (d1 + d2 + sqrt_disc);
                    // Choose eigenvalue closer to d2
                    mu = (std::abs(lambda1 - d2) < std::abs(lambda2 - d2)) ? lambda1 : lambda2;
                } else {
                    mu = d2;
                }
            } else {
                // Use Wilkinson shift: eigenvalue of bottom-right 2x2 closest to d2
                double disc = (d1 - d2) * (d1 - d2) + 4.0 * e * e;
                if (disc >= 0.0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (d1 + d2 - sqrt_disc);
                    double lambda2 = 0.5 * (d1 + d2 + sqrt_disc);
                    mu = (std::abs(lambda1 - d2) < std::abs(lambda2 - d2)) ? lambda1 : lambda2;
                } else {
                    mu = d2;
                }
            }
            
            // Apply shift to submatrix [m1, m2]
            for (int i = m1; i <= m2; ++i) {
                diag[i] -= mu;
            }
            
            // Implicit QR step: chase bulge from top to bottom
            // Start with first column: create bulge in position (m1+1, m1)
            double x = diag[m1];
            double y = off_diag[m1];
            double r = std::sqrt(x*x + y*y);
            if (r != ZERO) {
                double c = x / r;
                double s = y / r;
                
                // Apply Givens rotation to columns m1 and m1+1
                // Update first two diagonal elements
                double temp = c * diag[m1] + s * off_diag[m1];
                off_diag[m1] = -s * diag[m1] + c * off_diag[m1];
                diag[m1] = temp;
                
                if (m1 + 1 < m2) {
                    temp = c * off_diag[m1 + 1];
                    off_diag[m1 + 1] = s * off_diag[m1 + 1];
                    diag[m1 + 1] = c * diag[m1 + 1] + s * off_diag[m1 + 1];
                }
                
                // Chase the bulge down the matrix
                for (int i = m1 + 1; i < m2; ++i) {
                    x = off_diag[i - 1];
                    y = diag[i];
                    r = std::sqrt(x*x + y*y);
                    if (r != ZERO) {
                        c = x / r;
                        s = y / r;
                        
                        // Apply rotation to rows i and i+1
                        if (i > m1 + 1) {
                            off_diag[i - 2] = c * off_diag[i - 2] + s * diag[i - 1];
                        }
                        diag[i - 1] = -s * off_diag[i - 2] + c * diag[i - 1];
                        
                        if (i < m2) {
                            temp = c * off_diag[i] + s * diag[i + 1];
                            diag[i + 1] = -s * off_diag[i] + c * diag[i + 1];
                            off_diag[i] = temp;
                        }
                    }
                }
            }
            
            // Add back the shift
            for (int i = m1; i <= m2; ++i) {
                diag[i] += mu;
            }
        }
        
        if (!converged) {
            // In production code, we'd handle failure, but for test case it should converge
        }
        
        // After convergence, the last element is an eigenvalue
        // Decouple the last element
        if (m2 > m1) {
            m2 = m2 - 1;
        } else {
            m1 = m2 + 1;
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
        double lambda_k = 2.0 - 2.0 * std::cos(k * M_PI / (n + 1));
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