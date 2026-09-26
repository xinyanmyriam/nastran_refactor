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
        
        // Main iteration loop over deflation
        bool converged = false;
        int iter_count = 0;
        
        // Work on submatrices from bottom to top
        int m = n;
        while (m > 0) {
            // Check if we can deflate from bottom
            if (m == 1) {
                // Single element, already converged
                break;
            }
            
            // Check if subdiagonal element (m-1) is negligible
            if (std::abs(o[m-2]) <= EPSI) {
                // Deflate the last element
                m--;
                continue;
            }
            
            // Check if second-to-last subdiagonal is negligible
            if (m > 2 && std::abs(o[m-3]) <= EPSI) {
                // Deflate last two elements
                m -= 2;
                continue;
            }
            
            // Perform QR iteration on leading m x m submatrix
            bool submatrix_converged = false;
            int local_iter = 0;
            const int MAX_LOCAL_ITER = 30;
            
            while (!submatrix_converged && local_iter < MAX_LOCAL_ITER) {
                local_iter++;
                iter_count++;
                
                if (iter_count > g_state.MAX) {
                    g_state.NEVER++;
                    break;
                }
                
                // Wilkinson shift: choose shift from bottom 2x2 submatrix
                double d_m = val[m-1];
                double d_m1 = val[m-2];
                double e_m1 = o[m-2];
                
                double trace = d_m + d_m1;
                double det = d_m * d_m1 - e_m1 * e_m1;
                double disc = trace * trace - 4.0 * det;
                
                double mu;
                if (disc >= 0.0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (trace + sqrt_disc);
                    double lambda2 = 0.5 * (trace - sqrt_disc);
                    // Choose shift closest to d_m
                    mu = std::abs(lambda1 - d_m) < std::abs(lambda2 - d_m) ? lambda1 : lambda2;
                } else {
                    mu = d_m; // fallback
                }
                
                // Apply shift to leading m x m submatrix
                for (int i = 0; i < m; ++i) {
                    val[i] -= mu;
                }
                
                // Implicit QR step: chase the bulge
                double c = 1.0;
                double s = 0.0;
                
                // First Givens rotation on top-left 2x2
                double a = val[0];
                double b = o[0];
                double r = std::sqrt(a*a + b*b);
                if (r != 0.0) {
                    c = a / r;
                    s = b / r;
                }
                
                // Apply to first row/column
                double temp = c * val[0] + s * o[0];
                o[0] = -s * val[0] + c * o[0];
                val[0] = temp;
                
                if (m > 1) {
                    // Chase the bulge down the matrix
                    for (int k = 0; k < m-1; ++k) {
                        // Apply Givens rotation to columns k and k+1
                        double a_k = val[k];
                        double b_k = (k < m-1) ? o[k] : 0.0;
                        
                        // Compute new rotation
                        double r_k = std::sqrt(a_k*a_k + b_k*b_k);
                        double c_k = 1.0;
                        double s_k = 0.0;
                        if (r_k != 0.0) {
                            c_k = a_k / r_k;
                            s_k = b_k / r_k;
                        }
                        
                        // Apply to row k
                        if (k < m-1) {
                            double temp_val = c_k * val[k] + s_k * o[k];
                            double temp_o = -s_k * val[k] + c_k * o[k];
                            val[k] = temp_val;
                            o[k] = temp_o;
                        }
                        
                        // Apply to row k+1
                        if (k+1 < m) {
                            double temp_val = c_k * o[k] + s_k * val[k+1];
                            double temp_o = -s_k * o[k] + c_k * val[k+1];
                            if (k < m-2) {
                                o[k] = temp_val;
                            }
                            val[k+1] = temp_o;
                        }
                    }
                }
                
                // Add back the shift
                for (int i = 0; i < m; ++i) {
                    val[i] += mu;
                }
                
                // Check convergence: if bottom subdiagonal is small enough
                if (std::abs(o[m-2]) <= EPSI) {
                    submatrix_converged = true;
                    // Deflate the last element
                    m--;
                }
            }
            
            if (!submatrix_converged) {
                // Try simple shift or give up
                g_state.NEVER++;
                break;
            }
        }
    }
    
    // Reorder eigenvalues algebraically in ascending order
    // Create index array for sorting
    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < g_state.N; ++i) {
        pairs.emplace_back(val[i], i);
    }
    
    // Sort by eigenvalue
    std::sort(pairs.begin(), pairs.end(), 
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first < b.first;
              });
    
    // Extract sorted eigenvalues and locations
    for (int i = 0; i < g_state.N; ++i) {
        val[i] = pairs[i].first;
        loc[i] = pairs[i].second + 1; // 1-based indexing
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