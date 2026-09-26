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
    
    // Store eigenvalues as we deflate
    std::vector<double> eigenvals(n);
    int k = n - 1; // index where we'll store next eigenvalue
    
    // Main iteration loop
    int m1 = 0;  // current lower index (0-based)
    int m2 = n-1; // current upper index (0-based)
    
    for (int iter = 0; iter < max_iter && m1 <= m2; ) {
        // Check for convergence: if e[m2-1] is effectively zero, then d[m2] is an eigenvalue
        if (m2 == m1) {
            // Single element left
            eigenvals[k--] = d[m1];
            break;
        }
        
        // Check if bottom subdiagonal element is negligible
        if (std::abs(e[m2-1]) <= 1e-14 * (std::abs(d[m2-1]) + std::abs(d[m2]))) {
            // Deflate: d[m2] is an eigenvalue
            eigenvals[k--] = d[m2];
            m2--;
            continue;
        }
        
        // Check if top subdiagonal element is negligible
        if (m1 < m2-1 && std::abs(e[m1]) <= 1e-14 * (std::abs(d[m1]) + std::abs(d[m1+1]))) {
            // Deflate from top
            eigenvals[k--] = d[m1];
            m1++;
            continue;
        }
        
        // Compute Wilkinson shift (more stable than smallest diagonal element)
        // For 2x2 submatrix at bottom: [d[m2-1], e[m2-1]; e[m2-1], d[m2]]
        double a = d[m2-1];
        double b = e[m2-1];
        double c = d[m2];
        double discriminant = (a - c) * (a - c) + 4.0 * b * b;
        if (discriminant < 0.0) discriminant = 0.0;
        double sqrt_disc = std::sqrt(discriminant);
        double r1 = 0.5 * (a + c + sqrt_disc);
        double r2 = 0.5 * (a + c - sqrt_disc);
        double shift = (std::abs(r1 - c) < std::abs(r2 - c)) ? r1 : r2;
        
        // Apply shift to current submatrix
        for (int i = m1; i <= m2; i++) {
            d[i] -= shift;
        }
        
        // Perform implicit QR step (QL step for tridiagonal)
        // This is the standard bulge-chasing algorithm
        double f = d[m1];
        double g = e[m1];
        
        for (int i = m1; i < m2; i++) {
            // Compute Givens rotation parameters
            double r = std::sqrt(f * f + g * g);
            if (r == 0.0) {
                f = 0.0;
                g = 0.0;
                // Continue with next iteration
                if (i < m2 - 1) {
                    f = d[i+1];
                    g = e[i+1];
                }
                continue;
            }
            
            double c_val = f / r;
            double s_val = g / r;
            
            // Update d[i] and d[i+1]
            double temp = c_val * d[i] + s_val * e[i];
            e[i] = -s_val * d[i] + c_val * e[i];
            d[i] = temp;
            
            // Update e[i] (which becomes the new subdiagonal)
            if (i > m1) {
                e[i-1] = c_val * e[i-1] - s_val * d[i];
                d[i] = s_val * e[i-1] + c_val * d[i];
            }
            
            // Update next diagonal element
            if (i < m2 - 1) {
                f = d[i+1];
                g = e[i+1];
            }
        }
        
        // Handle last element
        if (m2 > m1) {
            double r = std::sqrt(d[m2] * d[m2] + e[m2-1] * e[m2-1]);
            if (r != 0.0) {
                double c_val = d[m2] / r;
                double s_val = e[m2-1] / r;
                d[m2] = r;
                e[m2-1] = 0.0;
            }
        }
        
        // Shift back
        if (shift != 0.0) {
            for (int i = m1; i <= m2; i++) {
                d[i] += shift;
            }
        }
        
        iter++;
    }
    
    // Collect remaining eigenvalues
    while (k >= 0 && m1 <= m2) {
        eigenvals[k--] = d[m1];
        m1++;
    }
    
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