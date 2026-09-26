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
        int iter = 0;
        int m = n; // current size of active submatrix
        
        while (m > 0 && iter < g_state.MAX) {
            // Deflation: find largest k such that o[k] is negligible
            int k = m - 1;
            while (k > 0 && std::abs(o[k-1]) > EPSI) {
                k--;
            }
            
            if (k == 0) {
                // No deflation possible, continue with full matrix
                k = 1;
            } else {
                // Deflate at position k
                m = k;
                iter = 0; // reset iteration count for smaller matrix
                continue;
            }
            
            // Check for convergence of bottom element
            if (m == 1) {
                break;
            }
            
            // Wilkinson shift: use bottom-right 2x2 submatrix
            double a = val[m-2];
            double b = o[m-2];
            double c = val[m-1];
            
            double p = 0.5 * (a - c);
            double r = std::sqrt(p * p + b * b);
            double s = c - b * b / (p + (p >= 0 ? r : -r));
            
            // Apply shift to submatrix [0..m-1]
            double shift = s;
            
            // Perform QR step with shift
            double x = val[0] - shift;
            double y = o[0];
            double z = 0.0;
            
            for (int i = 0; i < m - 1; ++i) {
                // Compute Givens rotation
                double r_val = std::sqrt(x * x + y * y);
                if (r_val == 0.0) {
                    r_val = EPSI;
                }
                double c_val = x / r_val;
                double s_val = y / r_val;
                
                // Apply rotation to diagonal element i
                double temp = c_val * val[i] + s_val * o[i];
                o[i] = -s_val * val[i] + c_val * o[i];
                val[i] = temp;
                
                // Apply rotation to next diagonal element
                if (i < m - 1) {
                    temp = c_val * o[i] + s_val * val[i+1];
                    val[i+1] = -s_val * o[i] + c_val * val[i+1];
                    o[i] = temp;
                }
                
                // Update for next iteration
                if (i < m - 2) {
                    x = val[i+1];
                    y = o[i+1];
                }
            }
            
            // Update for next iteration
            iter++;
            
            // Check convergence: if bottom off-diagonal is small, deflate
            if (m > 1 && std::abs(o[m-2]) <= EPSI * (std::abs(val[m-2]) + std::abs(val[m-1]))) {
                o[m-2] = 0.0;
                m--;
                iter = 0; // reset iteration count
            }
        }
        
        // Handle iteration failure
        if (iter >= g_state.MAX && m > 0) {
            g_state.NEVER++;
        }
    }
    
    // Reorder eigenvalues algebraically in ascending order
    // Simple bubble sort to match expected behavior
    for (int i = 0; i < g_state.N; ++i) {
        for (int j = 0; j < g_state.N - 1; ++j) {
            if (val[j] > val[j+1]) {
                std::swap(val[j], val[j+1]);
                std::swap(loc[j], loc[j+1]);
            }
        }
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