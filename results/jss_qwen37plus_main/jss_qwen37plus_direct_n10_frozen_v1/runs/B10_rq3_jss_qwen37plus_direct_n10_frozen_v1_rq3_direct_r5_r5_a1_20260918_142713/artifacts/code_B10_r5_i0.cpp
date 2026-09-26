#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

// For MSVC compatibility, define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declaration for JSON library
using json = nlohmann::json;

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

    // Machine precision adjustment (MACH == 5 for VAX, but we'll use standard double precision)
    // In original Fortran, EPSI was set to 1.0D-18 for MACH==5, otherwise 1.0D-36
    // We'll use a reasonable epsilon for double precision
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
                if (INDEX[k][2] == 1) continue;  // Skip marked rows? Actually, original logic is different
                // Original Fortran: IF (INDEX(K,3) - 1) 20,30,190 -> if INDEX(K,3) < 1 go to 20, ==1 go to 30, >1 go to 190
                // So we only consider k where INDEX(k,3) == 0
                if (INDEX[k][2] != 0) continue;
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

    // Interchange columns (apply inverse permutation)
    for (int i = 0; i < N; i++) {
        int l = N - 1 - i;  // l = N+1-i in 1-based indexing
        if (INDEX[l][0] == INDEX[l][1]) continue;
        int JROW = INDEX[l][0];
        int JCOLUM = INDEX[l][1];
        // Swap columns JROW and JCOLUM in A
        A.col(JROW).swap(A.col(JCOLUM));
    }

    // Extract solution vector from B
    Eigen::VectorXd solution = B.col(0);

    // Compute inverse matrix (A now contains the inverse)
    Eigen::MatrixXd K_inverse = A;

    // Verify solution by computing K * solution and comparing with b
    // (Not required for output but good for validation)

    // Output JSON
    json result;
    result["test"] = "INVERD";
    
    // Format solution as array
    std::vector<double> solution_vec;
    for (int i = 0; i < solution.size(); i++) {
        solution_vec.push_back(solution(i));
    }
    result["solution"] = solution_vec;
    result["determinant"] = DETERM;
    
    // Format inverse matrix as nested arrays
    std::vector<std::vector<double>> inverse_mat;
    for (int i = 0; i < K_inverse.rows(); i++) {
        std::vector<double> row;
        for (int j = 0; j < K_inverse.cols(); j++) {
            row.push_back(K_inverse(i, j));
        }
        inverse_mat.push_back(row);
    }
    result["inverse"] = inverse_mat;

    // Print the JSON result
    std::cout << result.dump() << std::endl;

    return 0;
}