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
    for (int iter = 0; iter < max_iter; ) {
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
        
        // Compute Wilkinson shift: eigenvalue of bottom-right 2x2 submatrix closest to d[m2]
        double a = d[m2-1];
        double b = e[m2-1];
        double c = d[m2];
        double discriminant = (a - c) * (a - c) + 4.0 * b * b;
        if (discriminant < 0.0) discriminant = 0.0;
        double sqrt_disc = std::sqrt(discriminant);
        double lambda1 = 0.5 * (a + c - sqrt_disc);
        double lambda2 = 0.5 * (a + c + sqrt_disc);
        double shift = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
        
        // Apply shift to current submatrix
        for (int i = m1; i <= m2; i++) {
            d[i] -= shift;
        }
        
        // Perform one implicit QL step (equivalent to QR for tridiagonal)
        // This is the standard bulge-chasing algorithm
        
        double f = d[m1];
        double g = e[m1];
        
        for (int i = m1; i < m2; i++) {
            // Compute Givens rotation parameters
            double r = std::sqrt(f * f + g * g);
            double c_val = f / r;
            double s = g / r;
            
            // Apply rotation to update d[i] and d[i+1]
            if (i > m1) {
                e[i-1] = r;
            }
            f = c_val * d[i] + s * e[i];
            e[i] = c_val * e[i] - s * d[i];
            d[i] = c_val * f + s * e[i];
            if (i < m2-1) {
                e[i] = c_val * e[i] - s * d[i+1];
                d[i+1] = s * f + c_val * d[i+1];
            }
        }
        
        e[m2-1] = f;
        d[m2] = f;
        
        // Shift back
        if (shift != 0.0) {
            for (int i = m1; i <= m2; i++) {
                d[i] += shift;
            }
        }
        
        // Check for convergence after the step
        bool converged = true;
        for (int i = m1; i < m2; i++) {
            if (std::abs(e[i]) > 1e-14 * (std::abs(d[i]) + std::abs(d[i+1]))) {
                converged = false;
                break;
            }
        }
        if (converged) {
            break;
        }
        
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