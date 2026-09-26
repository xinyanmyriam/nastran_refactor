#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Simulate the Fortran COMMON blocks and global parameters
struct GlobalState {
    int N = 0;           // order of the problem
    int MAX = 0;         // maximum iterations
    int NV = 0;          // number of eigenvectors
    int NE = 0;          // number of eigenvalues to print
    int NFOUND = 0;      // number of rigid modes found
    int IOPTN = 0;       // option flag
    double LFREQ = 0.0;  // low frequency bound
    double HFREQ = 0.0;  // high frequency bound
    double DLMDAS = 0.0; // shift parameter for MGIV method
    int NEVER = 0;       // iteration failure counter
};

GlobalState g_state;

// QR iteration for tridiagonal matrix eigenvalues
// Input: diagonal (val), off-diagonal (o) - both as vectors of doubles
// Output: eigenvalues in val (sorted ascending)
void QRITER(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr) {
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE = 1.0;
    
    // Ensure loc has correct size
    loc.resize(val.size());
    
    // Initialize loc array
    for (int i = 0; i < static_cast<int>(val.size()); ++i) {
        loc[i] = i + 1; // Fortran 1-based indexing
    }
    
    // Set up global state
    g_state.N = static_cast<int>(val.size());
    g_state.MAX = 100 * g_state.N;
    if (g_state.NV > g_state.N) g_state.NV = g_state.N;
    if (g_state.NE == 0) g_state.NE = g_state.N;
    if (g_state.NE < g_state.NV) g_state.NE = g_state.NV;
    
    // If QR == 0, perform actual QR iteration
    if (qr == 0) {
        int n = g_state.N;
        
        // Main iteration loop
        bool converged = false;
        int iter = 0;
        
        // Work arrays for Givens rotations
        std::vector<double> d = val;
        std::vector<double> e = o;
        
        while (iter < g_state.MAX && !converged) {
            iter++;
            
            // Find deflation points (small off-diagonal elements)
            int l = n - 1;
            while (l > 0 && std::abs(e[l-1]) > EPSI * (std::abs(d[l-1]) + std::abs(d[l]))) {
                l--;
            }
            
            if (l == 0) {
                // All off-diagonals are small enough
                converged = true;
                break;
            }
            
            // Check if bottom 2x2 block is decoupled
            if (l == n - 1) {
                // Try to deflate the last element
                if (std::abs(e[n-2]) <= EPSI * (std::abs(d[n-2]) + std::abs(d[n-1]))) {
                    e[n-2] = ZERO;
                    l = n - 2;
                }
            }
            
            // Wilkinson shift: eigenvalue of bottom 2x2 submatrix closest to d[n-1]
            double shift = ZERO;
            if (l < n - 1) {
                double a = d[n-2];
                double b = e[n-2];
                double c = d[n-1];
                double discriminant = (a - c) * (a - c) + 4.0 * b * b;
                if (discriminant >= 0.0) {
                    double sqrt_disc = std::sqrt(discriminant);
                    double lambda1 = 0.5 * (a + c - sqrt_disc);
                    double lambda2 = 0.5 * (a + c + sqrt_disc);
                    // Choose shift closest to d[n-1]
                    if (std::abs(lambda1 - d[n-1]) < std::abs(lambda2 - d[n-1])) {
                        shift = lambda1;
                    } else {
                        shift = lambda2;
                    }
                } else {
                    shift = d[n-1];
                }
            } else {
                shift = d[n-1];
            }
            
            // Apply shift to the matrix
            for (int i = 0; i < n; ++i) {
                d[i] -= shift;
            }
            
            // Implicit QR step with bulge chasing
            double f = d[0];
            double g = e[0];
            
            for (int k = 0; k < n - 1; ++k) {
                // Compute Givens rotation
                double r = std::sqrt(f * f + g * g);
                if (r == ZERO) {
                    r = EPSI;
                }
                double cs = f / r;
                double sn = g / r;
                
                // Apply to current row
                double temp = cs * d[k] + sn * e[k];
                e[k] = -sn * d[k] + cs * e[k];
                d[k] = temp;
                
                // Apply to next diagonal
                if (k < n - 2) {
                    temp = cs * e[k+1];
                    e[k+1] = -sn * e[k+1];
                    f = temp;
                } else {
                    f = cs * e[k];
                    e[k] = -sn * e[k];
                }
                
                // Apply to next off-diagonal
                if (k < n - 1) {
                    temp = cs * e[k] + sn * d[k+1];
                    d[k+1] = -sn * e[k] + cs * d[k+1];
                    e[k] = temp;
                }
                
                // Update g for next iteration
                if (k < n - 2) {
                    g = e[k+1];
                } else if (k == n - 2) {
                    g = e[k];
                }
            }
            
            // Add back the shift
            for (int i = 0; i < n; ++i) {
                d[i] += shift;
            }
            
            // Check convergence: if any off-diagonal is small enough, deflate
            bool deflated = false;
            for (int i = 0; i < n - 1; ++i) {
                if (std::abs(e[i]) <= EPSI * (std::abs(d[i]) + std::abs(d[i+1]))) {
                    e[i] = ZERO;
                    deflated = true;
                }
            }
            
            if (!deflated) {
                // No deflation occurred, continue
                continue;
            }
            
            // Count deflated eigenvalues
            int m = 0;
            for (int i = 0; i < n - 1; ++i) {
                if (e[i] == ZERO) {
                    m++;
                }
            }
            
            if (m == n - 1) {
                converged = true;
            }
        }
        
        // Copy results back
        val = d;
        o = e;
        
        // Handle iteration failure
        if (!converged) {
            g_state.NEVER++;
        }
    }
    
    // Reorder eigenvalues algebraically in ascending order
    // Create index array and sort by eigenvalues
    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < static_cast<int>(val.size()); ++i) {
        pairs.emplace_back(val[i], i + 1); // 1-based location
    }
    
    std::sort(pairs.begin(), pairs.end());
    
    // Extract sorted eigenvalues and locations
    for (int i = 0; i < static_cast<int>(val.size()); ++i) {
        val[i] = pairs[i].first;
        loc[i] = pairs[i].second;
    }
    
    // Handle rigid modes
    if (g_state.NFOUND > 0) {
        for (int i = 0; i < g_state.NFOUND; ++i) {
            val[i] = ZERO;
        }
    }
}

// Analytical solution for test case: lambda_k = 2 - 2*cos(k*pi/6), k=1..5
std::vector<double> analytical_eigenvalues() {
    const double PI = std::acos(-1.0);
    std::vector<double> result;
    for (int k = 1; k <= 5; ++k) {
        double lambda = 2.0 - 2.0 * std::cos(k * PI / 6.0);
        result.push_back(lambda);
    }
    return result;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    std::vector<double> val = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> o = {-1.0, -1.0, -1.0, -1.0}; // off-diagonals (subdiagonal only)
    std::vector<int> loc(5);
    
    // Set global parameters for test
    g_state.N = 5;
    g_state.NV = 0;
    g_state.NE = 0;
    g_state.NFOUND = 0;
    g_state.IOPTN = 0;
    g_state.LFREQ = 0.0;
    g_state.HFREQ = 0.0;
    g_state.DLMDAS = 0.0;
    
    // Perform QR iteration
    QRITER(val, o, loc, 0);
    
    // Sort eigenvalues (the algorithm should produce sorted, but ensure)
    std::sort(val.begin(), val.end());
    
    // Output JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::fixed << std::setprecision(12) << val[i];
        if (i < 4) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}