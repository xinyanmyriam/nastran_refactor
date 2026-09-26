#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Helper function to compute determinant of a matrix using LU decomposition
double computeDeterminant(const Eigen::MatrixXd& A) {
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    return lu.determinant();
}

// Helper function to compute inverse of a matrix
Eigen::MatrixXd computeInverse(const Eigen::MatrixXd& A) {
    return A.inverse();
}

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

    // Initialize variables as in Fortran
    int NDIM = 4;
    int N = 4;
    int M = 1; // one column of constants
    double DETERM = 1.0;
    int ISING = 1;
    
    // INDEX(N,3): index array for permutations
    // We'll use std::vector<std::array<int,3>> for INDEX
    std::vector<std::array<int,3>> INDEX(N);
    for (int j = 0; j < N; ++j) {
        INDEX[j][2] = 0; // INDEX(J,3) = 0
    }

    // Machine precision adjustment (MACH = 5 for VAX, but we'll use standard double)
    const double EPSI = 1.0e-18; // Adjusted for double precision

    // Main loop: I from 1 to N (0-indexed: 0 to N-1)
    for (int I = 0; I < N; ++I) {
        // SEARCH FOR PIVOT
        double AMAX = 0.0;
        int IROW = -1;
        int ICOLUM = -1;

        for (int J = 0; J < N; ++J) {
            if (INDEX[J][2] == 1) continue;
            for (int k = 0; k < N; ++k) {
                if (INDEX[k][2] == 1) continue;
                double abs_val = std::abs(A(J, k));
                if (abs_val > AMAX) {
                    AMAX = abs_val;
                    IROW = J;
                    ICOLUM = k;
                }
            }
        }

        if (IROW == -1 || ICOLUM == -1) {
            ISING = 2;
            break;
        }

        INDEX[ICOLUM][2]++; // INDEX(ICOLUM,3) = INDEX(ICOLUM,3) + 1
        INDEX[I][0] = IROW; // INDEX(I,1) = IROW
        INDEX[I][1] = ICOLUM; // INDEX(I,2) = ICOLUM

        // INTERCHANGE ROWS to put pivot element on diagonal
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

        // REDUCE NON-PIVOT ROWS
        for (int L1 = 0; L1 < N; ++L1) {
            if (L1 == ICOLUM) continue;
            double T = A(L1, ICOLUM);
            A(L1, ICOLUM) = 0.0;
            if (std::abs(T) < EPSI) continue;
            
            // Update row L1: A[L1,:] = A[L1,:] - T * A[ICOLUM,:]
            A.row(L1) = A.row(L1) - T * A.row(ICOLUM);
            if (M > 0) {
                B.row(L1) = B.row(L1) - T * B.row(ICOLUM);
            }
        }
    }

    // Check singularity
    for (int i = 0; i < N; ++i) {
        if (INDEX[i][2] != 1) {
            ISING = 2;
            break;
        }
    }
    if (ISING != 2) {
        ISING = 1;
    }

    // INTERCHANGE COLUMNS
    for (int i = 0; i < N; ++i) {
        int L = N - 1 - i; // L = N+1-I in Fortran (1-indexed)
        if (INDEX[L][0] == INDEX[L][1]) continue;
        int JROW = INDEX[L][0];
        int JCOLUM = INDEX[L][1];
        // Swap columns JROW and JCOLUM in A
        A.col(JROW).swap(A.col(JCOLUM));
    }

    // Extract solution vector from B
    Eigen::VectorXd solution = B.col(0);

    // Compute determinant using Eigen for verification
    double det_eigen = computeDeterminant(K);

    // Compute inverse using Eigen for verification
    Eigen::MatrixXd inv_eigen = computeInverse(K);

    // === OUTPUT IN PLAIN TEXT (replaces JSON) ===
    std::cout << "test: INVERD" << std::endl;

    // solution: [x0, x1, x2, x3]
    std::cout << "solution: [";
    for (int i = 0; i < solution.size(); ++i) {
        std::cout << std::setprecision(12) << solution(i);
        if (i < solution.size()-1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "determinant: " << std::setprecision(12) << det_eigen << std::endl;

    std::cout << "inverse:" << std::endl;
    for (int i = 0; i < inv_eigen.rows(); ++i) {
        for (int j = 0; j < inv_eigen.cols(); ++j) {
            std::cout << std::setprecision(12) << inv_eigen(i, j);
            if (j < inv_eigen.cols()-1) std::cout << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}