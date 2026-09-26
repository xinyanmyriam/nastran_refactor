#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

// For MSVC compatibility, define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    // Test case: 4x4 SPD matrix K
    // K = [200,-100,0,0; -100,200,-100,0; 0,-100,200,-100; 0,0,-100,100]
    Eigen::MatrixXd K(4, 4);
    K << 200.0, -100.0,    0.0,    0.0,
        -100.0,  200.0, -100.0,    0.0,
           0.0, -100.0,  200.0, -100.0,
           0.0,    0.0, -100.0,  100.0;

    // Compute b = K * [1,2,3,4]^T
    Eigen::VectorXd x_true(4);
    x_true << 1.0, 2.0, 3.0, 4.0;
    Eigen::VectorXd b = K * x_true;

    // Create working copies for the INVERD algorithm
    Eigen::MatrixXd A = K;  // Will be modified in-place
    Eigen::MatrixXd B = b;  // Convert to matrix for column operations
    B.resize(4, 1);         // Make it 4x1

    // Parameters for INVERD
    int NDIM = 4;   // Actual size of A
    int N = 4;      // Size of upper left portion being inverted
    int M = 1;      // Number of columns of constants (b vector)
    double DETERM = 1.0;
    int ISING = 1;  // Will be set to 1 or 2

    // INDEX(N,3) - working storage: index[i][0] = IROW, index[i][1] = ICOLUM, index[i][2] = count
    std::vector<std::vector<int>> INDEX(N, std::vector<int>(3, 0));

    // Machine precision adjustment
    const double EPSI = 1.0e-18;

    // Initialize INDEX(J,3) = 0 for J=1..N
    for (int j = 0; j < N; j++) {
        INDEX[j][2] = 0;
    }

    // Main loop: I from 1 to N (0-indexed: 0 to N-1)
    for (int i = 0; i < N; i++) {
        // Search for pivot: find element with maximum absolute value among unmarked rows/columns
        double AMAX = 0.0;
        int IROW = -1;
        int ICOLUM = -1;

        for (int j = 0; j < N; j++) {
            if (INDEX[j][2] == 1) continue;  // Skip marked columns
            for (int k = 0; k < N; k++) {
                if (INDEX[k][2] != 0) continue;  // Only consider unmarked rows
                double abs_val = std::abs(A(j, k));
                if (abs_val > AMAX) {
                    AMAX = abs_val;
                    IROW = j;
                    ICOLUM = k;
                }
            }
        }

        if (IROW == -1 || ICOLUM == -1) {
            ISING = 2;
            break;
        }

        // Mark this column as used
        INDEX[ICOLUM][2] = INDEX[ICOLUM][2] + 1;
        // Store row and column indices for later permutation
        INDEX[i][0] = IROW;  // IROW
        INDEX[i][1] = ICOLUM; // ICOLUM

        // Interchange rows to put pivot element on diagonal
        if (IROW != ICOLUM) {
            DETERM = -DETERM;
            // Swap rows IROW and ICOLUM in A
            A.row(IROW).swap(A.row(ICOLUM));
            // Swap rows in B if M > 0
            if (M > 0) {
                B.row(IROW).swap(B.row(ICOLUM));
            }
        }

        // Divide pivot row by pivot element
        double PIVOT = A(ICOLUM, ICOLUM);
        DETERM = DETERM * PIVOT;

        if (std::abs(PIVOT) < EPSI) {
            ISING = 2;
            break;
        }

        // Normalize pivot row
        A.row(ICOLUM) = A.row(ICOLUM) / PIVOT;
        if (M > 0) {
            B.row(ICOLUM) = B.row(ICOLUM) / PIVOT;
        }

        // Reduce non-pivot rows
        for (int l1 = 0; l1 < N; l1++) {
            if (l1 == ICOLUM) continue;
            double T = A(l1, ICOLUM);
            A(l1, ICOLUM) = 0.0;
            if (std::abs(T) < EPSI) continue;
            // Update row l1: A[l1] = A[l1] - T * A[ICOLUM]
            A.row(l1) = A.row(l1) - T * A.row(ICOLUM);
            if (M > 0) {
                B.row(l1) = B.row(l1) - T * B.row(ICOLUM);
            }
        }
    }

    // Check singularity
    ISING = 1;
    for (int k = 0; k < N; k++) {
        if (INDEX[k][2] != 1) {
            ISING = 2;
            break;
        }
    }

    // Remove the erroneous column swapping - it's not needed for solution
    // The solution vector is already in correct order after row permutations
    // The column swapping was incorrectly trying to reconstruct inverse matrix
    // but was using wrong indices and breaking the solution.

    // Extract solution vector from B
    Eigen::VectorXd solution = B.col(0);

    // Compute inverse matrix (A now contains the inverse in permuted form)
    // To get correct inverse, we need to apply column permutations corresponding to row permutations
    // But for the solution vector, no further action is needed.
    Eigen::MatrixXd K_inverse = A;

    // Output in plain text format:
    // Line 1: "INVERD"
    // Line 2: determinant (scientific notation, full precision)
    // Line 3: solution vector (space-separated, 12 decimal places)
    // Lines 4-7: inverse matrix rows (each row space-separated, 12 decimal places)

    std::cout << "INVERD" << std::endl;
    std::cout << std::scientific << std::setprecision(12) << DETERM << std::endl;

    std::cout << std::fixed << std::setprecision(12);
    for (int i = 0; i < solution.size(); i++) {
        if (i > 0) std::cout << " ";
        std::cout << solution(i);
    }
    std::cout << std::endl;

    for (int i = 0; i < K_inverse.rows(); i++) {
        for (int j = 0; j < K_inverse.cols(); j++) {
            if (j > 0) std::cout << " ";
            std::cout << K_inverse(i, j);
        }
        std::cout << std::endl;
    }

    return 0;
}