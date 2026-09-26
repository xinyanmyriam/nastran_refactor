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
        
        // Work arrays
        std::vector<double> d = val;
        std::vector<double> e = o;
        
        // Deflation indices
        int m = n;
        int l = 1;
        
        while (m >= l && iter < g_state.MAX) {
            // Look for deflation at bottom
            while (m > l) {
                double abs_e_m_minus_1 = std::abs(e[m-2]); // e[m-2] is off-diagonal between m-1 and m
                if (abs_e_m_minus_1 <= EPSI * (std::abs(d[m-1]) + std::abs(d[m-2]))) {
                    e[m-2] = ZERO;
                    break;
                }
                m--;
            }
            
            if (m == l) {
                // Single element
                break;
            }
            
            // Look for deflation at top
            while (l < m) {
                double abs_e_l_minus_1 = std::abs(e[l-1]); // e[l-1] is off-diagonal between l and l+1
                if (abs_e_l_minus_1 <= EPSI * (std::abs(d[l-1]) + std::abs(d[l]))) {
                    e[l-1] = ZERO;
                    l++;
                    break;
                }
                l++;
            }
            
            if (l == m) {
                // Single element
                break;
            }
            
            // Wilkinson shift: use bottom-right 2x2 submatrix
            double a = d[m-2];
            double b = e[m-2];
            double c = d[m-1];
            double shift = ZERO;
            
            if (m > l + 1) {
                // Wilkinson shift: eigenvalue of bottom-right 2x2 closest to c
                double p = 0.5 * (a - c);
                double r = std::sqrt(p * p + b * b);
                shift = c - b * b / (p + std::copysign(r, p));
            } else {
                // For 2x2, use simple shift
                shift = c;
            }
            
            // Apply shift to submatrix [l, m]
            for (int i = l; i <= m; ++i) {
                d[i-1] -= shift;
            }
            
            // Implicit QR step
            double f = d[l-1];
            double g = e[l-1];
            
            for (int k = l; k < m; ++k) {
                // Compute plane rotation
                double r = std::sqrt(f * f + g * g);
                if (r == ZERO) {
                    r = EPSI;
                }
                double cs = f / r;
                double sn = g / r;
                
                // Apply to diagonal
                if (k > l) {
                    e[k-2] = r;
                }
                f = cs * d[k-1] + sn * e[k-1];
                e[k-1] = cs * e[k-1] - sn * d[k-1];
                d[k-1] = cs * f + sn * e[k-1];
                g = -sn * f + cs * e[k-1];
                
                // Apply to next diagonal
                if (k < m-1) {
                    f = d[k];
                    d[k] = cs * f + sn * g;
                    g = -sn * f + cs * g;
                }
            }
            
            d[m-1] = f;
            e[m-2] = g;
            
            // Add back the shift
            for (int i = l; i <= m; ++i) {
                d[i-1] += shift;
            }
            
            iter++;
            
            // Check for convergence: if any off-diagonal is small enough, deflate
            bool deflated = false;
            for (int i = l; i < m; ++i) {
                if (std::abs(e[i-1]) <= EPSI * (std::abs(d[i-1]) + std::abs(d[i]))) {
                    e[i-1] = ZERO;
                    deflated = true;
                    break;
                }
            }
            
            if (!deflated) {
                // Continue with current submatrix
                continue;
            }
        }
        
        // Copy results back
        val = d;
        
        // Handle iteration failure
        if (iter >= g_state.MAX) {
            g_state.NEVER++;
        }
    }
    
    // Reorder eigenvalues algebraically in ascending order
    // Create index array for sorting
    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < g_state.N; ++i) {
        pairs.emplace_back(val[i], i + 1); // value, 1-based index
    }
    
    // Sort by eigenvalue
    std::sort(pairs.begin(), pairs.end(), 
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first < b.first;
              });
    
    // Extract sorted eigenvalues and locations
    for (int i = 0; i < g_state.N; ++i) {
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