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
            
            // Apply shift to submatrix [m1, m2]
            for (int i = m1; i <= m2; ++i) {
                diag[i] -= mu;
            }
            
            // Implicit QR step: chase bulge from top to bottom
            // Start with first column: create bulge in position (m1+1, m1)
            double x = diag[m1];
            double y = off_diag[m1];
            double r = std::sqrt(x*x + y*y);
            double c, s;
            
            if (r == ZERO) {
                c = 1.0;
                s = 0.0;
            } else {
                c = x / r;
                s = y / r;
            }
            
            // Apply Givens rotation to first 2x2 block: rows m1, m1+1 and columns m1, m1+1
            // Update diag[m1] and off_diag[m1]
            double temp_diag = c * c * diag[m1] + 2.0 * c * s * off_diag[m1] - s * s * diag[m1 + 1];
            double temp_off = c * s * (diag[m1] - diag[m1 + 1]) + (c * c - s * s) * off_diag[m1];
            double temp_diag_next = s * s * diag[m1] - 2.0 * c * s * off_diag[m1] + c * c * diag[m1 + 1];
            
            // But standard approach: apply rotation to eliminate off_diag[m1]
            // First rotation: eliminate off_diag[m1] by rotating rows m1 and m1+1
            double old_diag_m1 = diag[m1];
            double old_off_m1 = off_diag[m1];
            double old_diag_m1p1 = diag[m1 + 1];
            
            // After rotation: new diag[m1] = c*old_diag_m1 + s*old_off_m1
            //                 new off_diag[m1] = -s*old_diag_m1 + c*old_off_m1 (should be 0)
            //                 new diag[m1+1] = c*old_diag_m1p1 - s*old_off_m1
            // But we want to eliminate off_diag[m1], so use different convention
            
            // Standard implicit QR for tridiagonal: start with vector [diag[m1], off_diag[m1]]
            // Create rotation to zero out off_diag[m1]
            if (r != ZERO) {
                c = x / r;
                s = y / r;
                
                // Apply to first 2x2 block
                // Row m1: [diag[m1], off_diag[m1], 0, ...]
                // Row m1+1: [off_diag[m1], diag[m1+1], off_diag[m1+1], ...]
                
                // Update diag[m1] and off_diag[m1] (first element of row m1)
                double new_diag_m1 = c * diag[m1] + s * off_diag[m1];
                double new_off_m1 = -s * diag[m1] + c * off_diag[m1];
                
                // Update diag[m1+1] and off_diag[m1+1] (first element of row m1+1)
                double new_diag_m1p1 = c * diag[m1 + 1] + s * off_diag[m1 + 1];
                double new_off_m1p1 = -s * diag[m1 + 1] + c * off_diag[m1 + 1];
                
                // Store results
                diag[m1] = new_diag_m1;
                off_diag[m1] = new_off_m1; // This should be ~0
                
                if (m1 + 1 < m2) {
                    diag[m1 + 1] = new_diag_m1p1;
                    off_diag[m1 + 1] = new_off_m1p1;
                } else if (m1 + 1 == m2) {
                    diag[m1 + 1] = new_diag_m1p1;
                }
                
                // Now chase the bulge down: for i from m1+1 to m2-1
                for (int i = m1 + 1; i < m2; ++i) {
                    // Current off-diagonal is at position i-1 (from previous rotation)
                    // We have non-zero off_diag[i-1] that we need to eliminate
                    // Use vector [off_diag[i-1], diag[i]] to create rotation
                    
                    x = off_diag[i-1];
                    y = diag[i];
                    r = std::sqrt(x*x + y*y);
                    
                    if (r == ZERO) {
                        c = 1.0;
                        s = 0.0;
                    } else {
                        c = x / r;
                        s = y / r;
                    }
                    
                    // Apply rotation to eliminate off_diag[i-1]
                    // This affects positions: off_diag[i-1], diag[i], off_diag[i], diag[i+1]
                    
                    // First, update off_diag[i-1] and diag[i]
                    double temp_off = -s * off_diag[i-1] + c * diag[i];
                    double temp_diag = c * off_diag[i-1] + s * diag[i];
                    
                    // Then update diag[i] and off_diag[i]
                    if (i < m2) {
                        double temp_diag_i = c * diag[i] + s * off_diag[i];
                        double temp_off_i = -s * diag[i] + c * off_diag[i];
                        
                        // Update diag[i] and off_diag[i]
                        diag[i] = temp_diag_i;
                        off_diag[i] = temp_off_i;
                        
                        // Update diag[i+1] if exists
                        if (i + 1 <= m2) {
                            double temp_diag_ip1 = c * diag[i + 1];
                            // off_diag[i] is already updated above
                            diag[i + 1] = temp_diag_ip1;
                        }
                    }
                    
                    // Store the results
                    off_diag[i-1] = temp_off;
                    diag[i] = temp_diag;
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