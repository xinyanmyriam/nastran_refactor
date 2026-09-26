#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

// Portable M_PI definition (not standard C++)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Constants
const double EPSI = 1.0e-12;  // tighter tolerance
const double ZERO = 0.0;
const double ONE = 1.0;

// Global parameters (simplified for test)
int N = 5;                    // order of matrix
int MAX = 0;                   // max iterations (computed)
int NV = 0;                    // number of eigenvectors (not used)
int NE = 0;                    // number of eigenvalues to print (not used)
int NFOUND = 0;                // number of rigid modes (0 in test)
int IOPTN = 0;                 // option flag
double DLMDAS = 0.0;           // shift parameter
int NEVER = 0;                 // failure counter

// Safe access wrappers (defensive, though we'll fix indexing directly)
inline double get_val(const std::vector<double>& val, int i) {
    return (i >= 0 && i < static_cast<int>(val.size())) ? val[i] : ZERO;
}
inline double get_o(const std::vector<double>& o, int i) {
    return (i >= 0 && i < static_cast<int>(o.size())) ? o[i] : ZERO;
}

// QRITER implementation — fixed, safe, and convergent
void QRITER(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr) {
    // Initialize MAX = 100*N
    MAX = 100 * N;
    if (NV > N) NV = N;
    if (NE == 0) NE = N;
    if (NE < NV) NE = NV;

    NEVER = 0;

    // If qr != 0, just reorder (but in our test qr=0, so do full iteration)
    if (qr != 0) {
        goto reorder;
    }

    // Work on copy of val/o to avoid modifying originals during deflation
    std::vector<double> d = val;   // diagonal
    std::vector<double> e = o;     // subdiagonal (size N-1)

    // Deflation loop: work on active submatrix [low, high]
    int low = 0;
    int high = N - 1;

    while (high > low) {
        // Deflate small subdiagonal elements from bottom up
        bool deflated = false;
        for (int i = high - 1; i >= low; --i) {
            double abs_e_i = std::abs(e[i]);
            double threshold = EPSI * (std::abs(d[i]) + std::abs(d[i+1]));
            if (abs_e_i <= threshold) {
                e[i] = ZERO;
                deflated = true;
                break;
            }
        }
        if (!deflated) {
            // No deflation possible → run QR step on [low, high]
            int iter = 0;
            bool converged = false;

            while (iter < MAX && !converged) {
                // Check convergence: all subdiagonals in [low, high-1] must be ~0
                converged = true;
                for (int i = low; i < high; ++i) {
                    if (std::abs(e[i]) > EPSI * (std::abs(d[i]) + std::abs(d[i+1]))) {
                        converged = false;
                        break;
                    }
                }
                if (converged) break;

                // Wilkinson shift: eigenvalue of bottom 2x2 block
                double a = d[high-1];
                double b = e[high-1];
                double c = d[high];
                double shift = c;
                if (high > low + 1) {
                    double r = (a - c) / 2.0;
                    double s = r * r + b * b;
                    if (s > 0.0) {
                        shift = c - b * b / (r + (r > 0.0 ? 1.0 : -1.0) * std::sqrt(s));
                    }
                }

                // Apply shift to [low, high]
                for (int i = low; i <= high; ++i) {
                    d[i] -= shift;
                }

                // Bulge chase (single implicit QR step)
                double f = d[low];
                double g = e[low];
                for (int i = low; i < high; ++i) {
                    double r = std::sqrt(f*f + g*g);
                    if (r == ZERO) {
                        r = EPSI;
                    }
                    double cs = f / r;
                    double sn = g / r;

                    // Apply Givens rotation to rows i and i+1
                    if (i > low) {
                        e[i-1] = r;
                    }
                    f = cs * d[i] + sn * e[i];
                    e[i] = cs * e[i] - sn * d[i];
                    g = sn * d[i+1];
                    d[i+1] = cs * d[i+1];

                    if (i < high - 1) {
                        d[i] = f;
                        f = d[i+1];
                        g = e[i+1];
                    }
                }
                d[high] = f;

                // Unshift
                for (int i = low; i <= high; ++i) {
                    d[i] += shift;
                }

                iter++;
            }

            if (!converged) {
                NEVER++;
            }

            // Deflate converged bottom eigenvalue
            if (high > low && std::abs(e[high-1]) <= EPSI * (std::abs(d[high-1]) + std::abs(d[high]))) {
                e[high-1] = ZERO;
                high--;
            } else {
                // Failed to converge at bottom → try top
                if (std::abs(e[low]) <= EPSI * (std::abs(d[low]) + std::abs(d[low+1]))) {
                    e[low] = ZERO;
                    low++;
                } else {
                    // No deflation → break to avoid infinite loop
                    break;
                }
            }
        } else {
            // Find largest i < high where e[i] == 0 → split
            int i = high - 1;
            while (i > low && e[i] == ZERO) i--;
            if (i == low && e[i] == ZERO) {
                low++;
            }
            high = i;
        }
    }

    // Copy back converged diagonal
    val = d;

reorder:
    // Reorder eigenvalues algebraically in ascending order
    // Create index array [0,1,...,N-1]
    std::vector<int> idx(N);
    for (int i = 0; i < N; ++i) idx[i] = i;

    // Sort indices by val[i]
    std::sort(idx.begin(), idx.end(), [&](int i, int j) {
        return val[i] < val[j];
    });

    // Build sorted val and loc
    std::vector<double> sorted_val = val;
    for (int i = 0; i < N; ++i) {
        val[i] = sorted_val[idx[i]];
        loc[i] = idx[i] + 1; // 1-indexed position in original
    }

    // If rigid modes were found, set first NFOUND to zero (not applicable here)
    if (NFOUND > 0) {
        for (int i = 0; i < NFOUND; ++i) {
            val[i] = ZERO;
        }
    }
}

// Analytical solution for test case: lambda_k = 2 - 2*cos(k*pi/6), k=1..5
std::vector<double> analytical_eigenvalues() {
    std::vector<double> exact(5);
    for (int k = 1; k <= 5; ++k) {
        exact[k-1] = 2.0 - 2.0 * std::cos(k * M_PI / 6.0);
    }
    return exact;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    std::vector<double> val = {2.0, 2.0, 2.0, 2.0, 2.0}; // diagonal
    std::vector<double> o = {-1.0, -1.0, -1.0, -1.0};     // subdiagonal (size N-1 = 4)
    std::vector<int> loc(5);                               // location array

    // Set global parameters for test
    N = 5;
    NV = 0;
    NE = 0;
    NFOUND = 0;
    IOPTN = 0;
    DLMDAS = 0.0;

    // Call QRITER with QR=0 (full iteration)
    QRITER(val, o, loc, 0);

    // Output as JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setprecision(15) << val[i];
        if (i < 4) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}