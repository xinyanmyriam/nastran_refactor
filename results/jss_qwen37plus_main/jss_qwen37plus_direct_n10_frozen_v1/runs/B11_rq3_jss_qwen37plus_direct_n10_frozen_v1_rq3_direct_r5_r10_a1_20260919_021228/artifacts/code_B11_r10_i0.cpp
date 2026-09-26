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
    
    // Start with the full matrix [1..n]
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
        
        // Find shift: smallest diagonal element in current submatrix
        double shift = d[m1];
        for (int i = m1+1; i <= m2; i++) {
            if (std::abs(d[i]) < std::abs(shift)) {
                shift = d[i];
            }
        }
        
        // Apply shift to current submatrix
        for (int i = m1; i <= m2; i++) {
            d[i] -= shift;
        }
        
        // Perform one QR step (using implicit QL, which is equivalent to QR for tridiagonal)
        // This is the core Givens rotation based reduction
        
        double r = d[m1] * d[m1];
        double s = e[m1] / (r + e[m1]);
        double t = 0.0;
        double u = s * (d[m1] + d[m1+1]);
        d[m1] += u;
        
        if (m1 == m2-1) {
            // Only two elements
            d[m2] -= u;
            r = (1.0 - t) * e[m2-1];
            double ones = 1.0 - s;
            if (std::abs(ones) > 1e-10) {
                r = d[m2] * d[m2] / ones;
            }
            e[m2-1] = s * r;
            
            // Shift back
            if (shift != 0.0) {
                d[m1] += shift;
                d[m2] += shift;
            }
            break;
        }
        
        // More than two elements: iterate through the submatrix
        int m1p1 = m1 + 1;
        for (int i = m1p1; i < m2; i++) {
            double g = d[i] - u;
            r = (1.0 - t) * e[i-1];
            double ones = 1.0 - s;
            if (std::abs(ones) > 1e-10) {
                r = g * g / ones;
            }
            r += e[i];
            e[i-1] = s * r;
            
            if (std::abs(e[i-1]) <= 1e-14 * (std::abs(d[i-1]) + std::abs(d[i]))) {
                // Deflation point found
                m1 = i;
                break;
            }
            
            t = s;
            s = e[i] / r;
            u = s * (g + d[i+1]);
            d[i] = u + g;
        }
        
        // Handle last element
        if (m1 < m2) {
            d[m2] -= u;
            r = (1.0 - t) * e[m2-1];
            double ones = 1.0 - s;
            if (std::abs(ones) > 1e-10) {
                r = d[m2] * d[m2] / ones;
            }
            e[m2-1] = s * r;
        }
        
        // Shift back
        if (shift != 0.0) {
            for (int i = m1; i <= m2; i++) {
                d[i] += shift;
            }
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