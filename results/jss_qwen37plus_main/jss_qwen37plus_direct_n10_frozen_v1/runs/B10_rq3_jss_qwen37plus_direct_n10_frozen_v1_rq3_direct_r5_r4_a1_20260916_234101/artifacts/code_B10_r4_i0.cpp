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
    Eigen::MatrixXd A = K;
    Eigen::MatrixXd B = b;
    B.conservativeResize(4, 1); // Ensure B is 4x1

    int NDIM = 4;
    int N = 4;
    int M = 1;
    double DETERM = 0.0;
    int ISING = 0;
    
    // INDEX(N,3) - we'll use a 2D vector: index[i][j] where i=0..N-1, j=0..2
    std::vector<std::vector<int>> INDEX(N, std::vector<int>(3, 0));

    // Machine precision adjustment (MACH = 5 for VAX in original, but we'll use standard double precision)
    // EPSI = 1.0e-18 for double precision
    const double EPSI = 1.0e-18;

    // Initialize DETERM and INDEX
    DETERM = 1.0;
    ISING = 0;
    for (int j = 0; j < N; j++) {
        INDEX[j][2] = 0;
    }

    // Main loop: I from 0 to N-1 (Fortran 1-based -> C++ 0-based)
    for (int I = 0; I < N; I++) {
        // SEARCH FOR PIVOT
        double AMAX = 0.0;
        int IROW = -1;
        int ICOLUM = -1;

        for (int J = 0; J < N; J++) {
            if (INDEX[J][2] == 1) continue;
            for (int K_idx = 0; K_idx < N; K_idx++) {
                if (INDEX[K_idx][2] != 0) continue; // Fortran: IF (INDEX(K,3) - 1) 20,30,190 -> jump to 190 if >1
                if (std::abs(A(J, K_idx)) <= AMAX) continue;
                IROW = J;
                ICOLUM = K_idx;
                AMAX = std::abs(A(J, K_idx));
            }
        }

        if (IROW == -1 || ICOLUM == -1) {
            // Singular matrix
            ISING = 2;
            break;
        }

        INDEX[ICOLUM][2] = INDEX[ICOLUM][2] + 1;
        INDEX[I][0] = IROW;  // INDEX(I,1) = IROW
        INDEX[I][1] = ICOLUM; // INDEX(I,2) = ICOLUM

        // INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
        if (IROW != ICOLUM) {
            DETERM = -DETERM;
            // Swap rows IROW and ICOLUM in A
            A.row(IROW).swap(A.row(ICOLUM));
            // Swap rows IROW and ICOLUM in B (if M > 0)
            if (M > 0) {
                B.row(IROW).swap(B.row(ICOLUM));
            }
        }

        // DIVIDE PIVOT ROW BY PIVOT ELEMENT
        double PIVOT = A(ICOLUM, ICOLUM);
        DETERM = DETERM * PIVOT;

        if (std::abs(PIVOT) < EPSI) {
            ISING = 2;
            break;
        }

        A(ICOLUM, ICOLUM) = 1.0;
        // Divide entire row ICOLUM by PIVOT
        A.row(ICOLUM) = A.row(ICOLUM) / PIVOT;
        if (M > 0) {
            B.row(ICOLUM) = B.row(ICOLUM) / PIVOT;
        }

        // REDUCE NON PIVOT ROWS
        for (int L1 = 0; L1 < N; L1++) {
            if (L1 == ICOLUM) continue;
            double T = A(L1, ICOLUM);
            A(L1, ICOLUM) = 0.0;
            if (std::abs(T) < EPSI) continue;
            
            // Update row L1: A[L1] = A[L1] - T * A[ICOLUM]
            A.row(L1) = A.row(L1) - T * A.row(ICOLUM);
            if (M > 0) {
                B.row(L1) = B.row(L1) - T * B.row(ICOLUM);
            }
        }
    }

    // Check singularity
    if (ISING == 0) {
        ISING = 1;
        for (int K = 0; K < N; K++) {
            if (INDEX[K][2] != 1) {
                ISING = 2;
                break;
            }
        }
    }

    // INTERCHANGE COLUMNS
    for (int I = 0; I < N; I++) {
        int L = N - 1 - I; // Fortran: L = N + 1 - I
        if (INDEX[L][0] == INDEX[L][1]) continue;
        int JROW = INDEX[L][0];
        int JCOLUM = INDEX[L][1];
        // Swap columns JROW and JCOLUM in A
        A.col(JROW).swap(A.col(JCOLUM));
    }

    // Prepare output
    std::vector<double> solution;
    for (int i = 0; i < 4; i++) {
        solution.push_back(B(i, 0));
    }

    // Compute determinant using Eigen's built-in method for verification
    double det_eigen = K.determinant();

    // Compute inverse using Eigen for verification
    Eigen::MatrixXd K_inv = K.inverse();

    // Build JSON output
    json result;
    result["test"] = "INVERD";
    result["solution"] = solution;
    result["determinant"] = DETERM;
    result["inverse"] = json::array();
    
    for (int i = 0; i < 4; i++) {
        json row = json::array();
        for (int j = 0; j < 4; j++) {
            row.push_back(A(i, j));
        }
        result["inverse"].push_back(row);
    }

    // Print the JSON
    std::cout << result.dump(2) << std::endl;

    return 0;
}