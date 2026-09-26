#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Simplified version of the Fortran QRITER subroutine for tridiagonal matrices
// We focus only on the core QR iteration logic for eigenvalues of a symmetric tridiagonal matrix
// The test case is 5x5 with diagonal = [2,2,2,2,2] and off-diagonal = [-1,-1,-1,-1]

// Helper function to compute eigenvalues of a symmetric tridiagonal matrix using QR iteration
std::vector<double> qriter_eigenvalues(
    std::vector<double> diag,      // diagonal elements (size n)
    std::vector<double> subdiag    // subdiagonal elements (size n-1), same as superdiagonal
) {
    const int n = static_cast<int>(diag.size());
    if (n == 0) return {};
    if (n == 1) return {diag[0]};
    
    // Work arrays - we'll modify diag and subdiag in-place
    std::vector<double> d = diag;
    std::vector<double> e = subdiag; // e[i] is the subdiagonal element between rows i and i+1
    
    // Maximum iterations: 100 * n as in Fortran
    const int max_iter = 100 * n;
    
    // Start with the full matrix [0..n-1]
    int m1 = 0;  // current lower index (0-based)
    int m2 = n-1; // current upper index (0-based)
    
    // Main iteration loop
    for (int iter = 0; iter < max_iter && m1 <= m2; ) {
        // Check for convergence: if e[m2-1] is effectively zero, then d[m2] is an eigenvalue
        if (m2 == m1) {
            // Single element left
            break;
        }
        
        // Check if bottom subdiagonal element is negligible
        if (std::abs(e[m2-1]) <= 1e-14 * (std::abs(d[m2-1]) + std::abs(d[m2]))) {
            // Deflate: d[m2] is an eigenvalue
            m2--;
            continue;
        }
        
        // Check if top subdiagonal element is negligible
        if (m1 < m2-1 && std::abs(e[m1]) <= 1e-14 * (std::abs(d[m1]) + std::abs(d[m1+1]))) {
            // Deflate from top
            m1++;
            continue;
        }
        
        // Compute Wilkinson shift from bottom-right 2x2 submatrix
        double a = d[m2-1];
        double b = e[m2-1];
        double c = d[m2];
        double f = (a - c) / 2.0;
        double g = std::sqrt(f*f + b*b);
        double shift = c - b*b/(f + std::copysign(g, f));
        
        // Apply shift to current submatrix
        for (int i = m1; i <= m2; i++) {
            d[i] -= shift;
        }
        
        // Perform implicit QL step (equivalent to QR for symmetric matrices)
        // This is the standard Givens rotation based reduction
        
        // Initialize the first rotation
        double c_prev = 1.0;
        double s_prev = 0.0;
        double d_temp = d[m1];
        
        // Sweep from top to bottom
        for (int i = m1; i < m2; i++) {
            // Compute Givens rotation to eliminate e[i]
            double f = c_prev * e[i] + s_prev * d[i+1];
            double g = s_prev * e[i] - c_prev * d[i+1];
            
            // Update d[i+1]
            d[i+1] = c_prev * d[i+1] - s_prev * e[i];
            
            // Compute new rotation
            double r = std::sqrt(f*f + g*g);
            if (r == 0.0) {
                c_prev = 1.0;
                s_prev = 0.0;
            } else {
                c_prev = f / r;
                s_prev = g / r;
            }
            
            // Update e[i] and d[i]
            if (i > m1) {
                e[i-1] = r;
            }
            d[i] = r;
            
            // Apply rotation to next element
            if (i < m2-1) {
                double temp = c_prev * d[i+1] + s_prev * e[i+1];
                e[i+1] = -s_prev * d[i+1] + c_prev * e[i+1];
                d[i+1] = temp;
            }
        }
        
        // Handle the last element
        if (m2 > m1) {
            e[m2-1] = s_prev * d[m2];
            d[m2] = c_prev * d[m2];
        }
        
        // Shift back
        if (shift != 0.0) {
            for (int i = m1; i <= m2; i++) {
                d[i] += shift;
            }
        }
        
        // Check for convergence after the step
        iter++;
    }
    
    // Collect all eigenvalues (they should be in d[0..n-1] now)
    std::vector<double> eigenvals = d;
    
    // Sort in ascending order
    std::sort(eigenvals.begin(), eigenvals.end());
    
    return eigenvals;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2]
    // off-diagonal = [-1,-1,-1,-1] (subdiagonal and superdiagonal)
    
    std::vector<double> diag = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> subdiag = {-1.0, -1.0, -1.0, -1.0};
    
    // Compute eigenvalues using our QR iteration
    std::vector<double> eigenvals = qriter_eigenvalues(diag, subdiag);
    
    // Verify with analytical solution: lambda_k = 2 - 2*cos(k*pi/6), k=1..5
    // But we just output the computed values as required
    
    // Format output as JSON
    std::ostringstream json;
    json << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < static_cast<int>(eigenvals.size()); ++i) {
        if (i > 0) json << ",";
        json << std::fixed << std::setprecision(12) << eigenvals[i];
    }
    json << "]}";
    
    std::cout << json.str() << std::endl;
    
    return 0;
}