#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Eigen;

// ---- Wedge node coordinates (test case) ----
// N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
static double NODES[6][3] = {
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}
};

// Tetrahedron mapping for wedge (1-indexed node numbers)
static int MTET[12][4] = {
    {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
    {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
    {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
};

// Material
static double E_mat = 200e9;
static double NU = 0.3;

// Build 6x6 material matrix G (isotropic)
MatrixXd buildG() {
    MatrixXd G = MatrixXd::Zero(6,6);
    double temp1 = (1.0+NU)*(1.0-2.0*NU);
    double c1 = E_mat*(1.0-NU)/temp1;
    double c2 = E_mat*NU/temp1;
    double GG = E_mat/(2.0*(1.0+NU));
    G(0,0)=c1; G(1,1)=c1; G(2,2)=c1;
    G(0,1)=c2; G(0,2)=c2;
    G(1,0)=c2; G(1,2)=c2;
    G(2,0)=c2; G(2,1)=c2;
    G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;
    return G;
}

// Compute 6x3 C matrix for a tetra node given the inverse H matrix
// H is 4x4: rows are [1 x y z] for each of 4 nodes. Hinv = inverse.
// The shape function derivatives: dNi/dx = Hinv(1,i), dNi/dy = Hinv(2,i), dNi/dz = Hinv(3,i)
// C_i (6x3) built from these derivatives.
MatrixXd buildCi(const MatrixXd& Hinv, int i) {
    // i is 0-based node index within tet
    double dx = Hinv(1,i);
    double dy = Hinv(2,i);
    double dz = Hinv(3,i);
    MatrixXd C = MatrixXd::Zero(6,3);
    C(0,0)=dx;
    C(1,1)=dy;
    C(2,2)=dz;
    C(3,0)=dy; C(3,1)=dx;
    C(4,1)=dz; C(4,2)=dy;
    C(5,0)=dz; C(5,2)=dx;
    return C;
}

int main() {
    MatrixXd G = buildG();
    MatrixXd Ke = MatrixXd::Zero(18,18);

    // For each tetrahedron
    for (int t = 0; t < 12; t++) {
        // Gather 4 nodes
        MatrixXd H = MatrixXd::Zero(4,4);
        for (int i = 0; i < 4; i++) {
            int nd = MTET[t][i]-1;
            H(i,0)=1.0;
            H(i,1)=NODES[nd][0];
            H(i,2)=NODES[nd][1];
            H(i,3)=NODES[nd][2];
        }
        double detH = H.determinant();
        MatrixXd Hinv = H.inverse();
        double HDETER = std::abs(detH);

        // Build C matrices (6x3 each) for the 4 tet nodes
        MatrixXd C[4];
        for (int i = 0; i < 4; i++) C[i] = buildCi(Hinv, i);

        // Weight: for wedge, HDETER/36, and if IOPT<=16 (t=0..5) multiply by 2
        double w = HDETER/36.0;
        if (t < 6) w *= 2.0;

        // Tet stiffness 12x12: K_ij = w * C_i^T G C_j
        // Expand to 18x18 using wedge node indices
        for (int i = 0; i < 4; i++) {
            int ni = MTET[t][i]-1;
            for (int j = 0; j < 4; j++) {
                int nj = MTET[t][j]-1;
                MatrixXd Kij = w * (C[i].transpose() * G * C[j]);
                Ke.block<3,3>(3*ni, 3*nj) += Kij;
            }
        }
    }

    // Print as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 18; i++) {
        std::cout << "[";
        for (int j = 0; j < 18; j++) {
            std::cout << std::scientific << std::setprecision(6) << Ke(i,j);
            if (j < 17) std::cout << ",";
        }
        std::cout << "]";
        if (i < 17) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}