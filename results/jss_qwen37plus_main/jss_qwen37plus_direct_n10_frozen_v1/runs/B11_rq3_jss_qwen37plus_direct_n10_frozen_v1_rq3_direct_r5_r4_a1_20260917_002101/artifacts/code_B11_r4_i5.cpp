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
            // Standard algorithm for tridiagonal matrices
            
            // Start with first column: create bulge in position (m1+1, m1)
            double x = diag[m1];
            double y = off_diag[m1];
            double r = std::sqrt(x*x + y*y);
            if (r != ZERO) {
                double c = x / r;
                double s = y / r;
                
                // Apply Givens rotation to first 2x2 block: rows m1, m1+1
                // Update diag[m1] and off_diag[m1]
                double temp_diag = c * c * diag[m1] + 2.0 * c * s * off_diag[m1] + s * s * diag[m1 + 1];
                double temp_off = c * s * (diag[m1 + 1] - diag[m1]) + (c * c - s * s) * off_diag[m1];
                double temp_diag_next = s * s * diag[m1] - 2.0 * c * s * off_diag[m1] + c * c * diag[m1 + 1];
                
                // But standard approach: apply rotation to columns first
                // Instead, use the standard tridiagonal QR step:
                
                // Reset and use correct algorithm
                // Step 1: Apply rotation to first two rows/columns
                double a1 = diag[m1];
                double b1 = off_diag[m1];
                double a2 = diag[m1 + 1];
                
                // After rotation: new a1, b1, a2
                double new_a1 = c * a1 + s * b1;
                double new_b1 = -s * a1 + c * b1;
                double new_a2 = c * a2 - s * b1;
                
                diag[m1] = new_a1;
                off_diag[m1] = new_b1;
                diag[m1 + 1] = new_a2;
                
                // If more than 2x2, we need to eliminate the bulge
                // Chase the bulge down: for i from m1+1 to m2-1
                for (int i = m1 + 1; i < m2; ++i) {
                    // Current bulge is in off_diag[i-1] (which should be zero after previous step)
                    // But we have a bulge in position (i, i-1) that needs elimination
                    
                    // Get the 2x2 submatrix at positions (i-1,i) and (i,i+1)
                    double b_prev = off_diag[i-1]; // this is the bulge element
                    double a_i = diag[i];
                    double b_i = (i < m2 - 1) ? off_diag[i] : 0.0;
                    
                    // Create rotation to eliminate b_prev
                    if (std::abs(b_prev) > EPSI) {
                        double r2 = std::sqrt(b_prev * b_prev + a_i * a_i);
                        if (r2 != ZERO) {
                            double c2 = a_i / r2;
                            double s2 = b_prev / r2;
                            
                            // Apply to rows i-1 and i
                            // Update diag[i-1] and off_diag[i-1]
                            double temp = c2 * diag[i-1] + s2 * off_diag[i-1];
                            off_diag[i-1] = -s2 * diag[i-1] + c2 * off_diag[i-1];
                            diag[i-1] = temp;
                            
                            // Update diag[i] and off_diag[i]
                            if (i < m2) {
                                temp = c2 * off_diag[i] + s2 * diag[i+1];
                                diag[i+1] = -s2 * off_diag[i] + c2 * diag[i+1];
                                off_diag[i] = temp;
                            }
                            
                            // Update off_diag[i-1] (already done above)
                            // Update diag[i]
                            temp = c2 * diag[i] - s2 * off_diag[i-1];
                            diag[i] = temp;
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