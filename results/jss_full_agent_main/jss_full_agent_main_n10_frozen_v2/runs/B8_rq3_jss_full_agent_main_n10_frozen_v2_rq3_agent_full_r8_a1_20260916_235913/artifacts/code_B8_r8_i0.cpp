#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Eigen;

// Global-ish state mimicking Fortran COMMON
struct State {
    double E = 200e9, nu = 0.3, GG;
    double HDETER;
    int NPVT = 0;
    double K[18][18]; // global 18x18 stiffness
};

// SAXB: cross product R = R12 x R13
void SAXB(const double* a, const double* b, double* r) {
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];
}

double SADOTB(const double* a, const double* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Invert 4x4 matrix H (row-major, 4x4), return determinant
// Fortran INVERD(4,H,4,DUM,0,HDETER,ISING,TEMP)
// We'll just use Eigen inverse and determinant.
double invert4(const double H[16], double Hinv[16], int& ising) {
    Matrix4d m;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            m(i,j) = H[i*4+j];
    double det = m.determinant();
    ising = 0;
    if (std::abs(det) < 1e-300) { ising = 2; return det; }
    Matrix4d inv = m.inverse();
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            Hinv[i*4+j] = inv(i,j);
    return det;
}

// GMMATD: matrix multiply. Fortran GMMATD(A, nrowA, ncolA, transA, B, nrowB, ncolB, transB, C)
// transA=1 means use A^T, transB=1 means use B^T
// A is stored column-major in Fortran. We'll implement with explicit dims.
// Actually GMMATD signature: GMMATD(A, L, M, TRANA, B, N, K, TRANB, C)
// computes C = op(A) * op(B) where op(A) is L x M, op(B) is N x K... 
// Let me use the standard NASTRAN GMMATD: C(L,K) = A(L,M) * B(M,K) with transpose flags.
// The call GMMATD(T(1),3,3,1, C(NPOINT),6,3,1, CT(1)) means:
//   A = T, dims 3x3, transA=1 -> A^T (3x3)
//   B = C, dims 6x3, transB=1 -> B^T (3x6)
//   result CT = A^T * B^T = 3x6
// So GMMATD(A, L, M, trA, B, N, K, trB, C): op(A) is LxM, op(B) is NxK, result LxK.
// When trA=0: A is LxM. When trA=1: A is MxL, op(A)=A^T is LxM.
// When trB=0: B is NxK. When trB=1: B is KxN, op(B)=B^T is NxK.
// Result C is LxK.

void GMMATD(const double* A, int L, int M, int trA,
            const double* B, int N, int K, int trB,
            double* C) {
    // op(A): L x M
    // op(B): N x K  (N should equal M)
    // C: L x K
    // Fortran column-major storage
    auto getA = [&](int i, int j) -> double {
        // op(A)(i,j)
        if (trA == 0) return A[j*L + i];      // A is LxM col-major
        else          return A[i*M + j];      // A is MxL col-major, A^T(i,j)=A(j,i)
    };
    auto getB = [&](int i, int j) -> double {
        // op(B)(i,j)
        if (trB == 0) return B[j*N + i];      // B is NxK col-major
        else          return B[i*K + j];      // B is KxN col-major, B^T(i,j)=B(j,i)
    };
    for (int i = 0; i < L; i++)
        for (int j = 0; j < K; j++) {
            double s = 0;
            for (int p = 0; p < M; p++)
                s += getA(i,p) * getB(p,j);
            C[j*L + i] = s; // col-major
        }
}

// TRANSD: get transformation matrix T (3x3) for coordinate system CSID
// For simplicity, assume identity (CSID=0 in test case)
void TRANSD(int csid, double T[9]) {
    // identity
    for (int i = 0; i < 9; i++) T[i] = 0;
    T[0] = T[4] = T[8] = 1.0;
}

int main() {
    // Test case: wedge with 6 nodes
    // N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    double E = 200e9, nu = 0.3;
    double GG = E / (2.0*(1.0+nu));
    
    // ECPT for wedge (1-indexed in Fortran, we use 0-indexed)
    // ECPT(1)=EL ID, ECPT(2)=MAT-ID, ECPT(3..8)=GRID 1..6
    // ECPT(9)=CSID-1, ECPT(10..12)=X1,Y1,Z1
    // ECPT(13)=CSID-2, ECPT(14..16)=X2,Y2,Z2
    // ECPT(17)=CSID-3, ECPT(18..20)=X3,Y3,Z3
    // ECPT(21)=CSID-4, ECPT(22..24)=X4,Y4,Z4
    // ECPT(25)=CSID-5, ECPT(26..28)=X5,Y5,Z5
    // ECPT(29)=CSID-6, ECPT(30..32)=X6,Y6,Z6
    // ECPT(33)=ELTEMP
    double ECPT[100];
    for (int i = 0; i < 100; i++) ECPT[i] = 0;
    // 0-indexed: ECPT[0]=EL ID, ECPT[1]=MAT-ID
    ECPT[0] = 1;   // EL ID
    ECPT[1] = 1;   // MAT-ID
    ECPT[2] = 1;   // GRID-1
    ECPT[3] = 2;   // GRID-2
    ECPT[4] = 3;   // GRID-3
    ECPT[5] = 4;   // GRID-4
    ECPT[6] = 5;   // GRID-5
    ECPT[7] = 6;   // GRID-6
    // CSID-1 = 0, X1,Y1,Z1
    ECPT[8] = 0;   // CSID-1
    ECPT[9] = 0; ECPT[10] = 0; ECPT[11] = 0;  // N1
    ECPT[12] = 0;  // CSID-2
    ECPT[13] = 1; ECPT[14] = 0; ECPT[15] = 0; // N2
    ECPT[16] = 0;  // CSID-3
    ECPT[17] = 0; ECPT[18] = 1; ECPT[19] = 0; // N3
    ECPT[20] = 0;  // CSID-4
    ECPT[21] = 0; ECPT[22] = 0; ECPT[23] = 1; // N4
    ECPT[24] = 0;  // CSID-5
    ECPT[25] = 1; ECPT[26] = 0; ECPT[27] = 1; // N5
    ECPT[28] = 0;  // CSID-6
    ECPT[29] = 0; ECPT[30] = 1; ECPT[31] = 1; // N6
    ECPT[32] = 0;  // ELTEMP
    
    // M matrix for wedge (12 tetrahedra), 1-indexed in Fortran
    int M[12][4] = {
        {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
        {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
        {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
    };
    
    // Global stiffness 18x18
    double K[18][18];
    for (int i = 0; i < 18; i++)
        for (int j = 0; j < 18; j++)
            K[i][j] = 0;
    
    // For each tetrahedron
    for (int itet = 0; itet < 12; itet++) {
        int IOPT = itet + 1 + 10; // IOPT = I + 10, I is 1-indexed -> itet+1+10
        // Actually IOPT = I + 10 where I is the loop index (1..12)
        // So IOPT ranges 11..22
        
        // Build tetrahedron ECPT
        // NECPT(J+2) = NECPT(KPOINT+52) for J=1..4 -> grid SILs
        // KPOINT = 4*KPOINT + NGRIDS - 3, NGRIDS=6
        // JPOINT = 4*J + 2
        // NECPT(JPOINT+1..4) = NECPT(KPOINT+52..55) -> coords
        
        // We need to extract the 4 nodes of this tetrahedron
        int nodes[4];
        double coords[4][3];
        for (int j = 0; j < 4; j++) {
            int kpoint = M[itet][j]; // 1-indexed grid
            nodes[j] = kpoint;
            // coords: ECPT for grid kpoint
            // grid 1: ECPT[9..11], grid 2: ECPT[13..15], etc.
            // grid g: ECPT[9 + (g-1)*4 .. 11 + (g-1)*4]
            int base = 9 + (kpoint-1)*4;
            coords[j][0] = ECPT[base];
            coords[j][1] = ECPT[base+1];
            coords[j][2] = ECPT[base+2];
        }
        
        // Build H matrix (4x4)
        double H[16];
        for (int j = 0; j < 4; j++) {
            H[j*4+0] = 1.0;
            H[j*4+1] = coords[j][0];
            H[j*4+2] = coords[j][1];
            H[j*4+3] = coords[j][2];
        }
        
        // Invert H
        double Hinv[16];
        int ising = 0;
        double HDETER = invert4(H, Hinv, ising);
        if (ising == 2) continue;
        
        HDETER = std::abs(HDETER);
        
        // Build G matrix (6x6)
        double G[36];
        for (int i = 0; i < 36; i++) G[i] = 0;
        double TEMP1 = (1.0+nu)*(1.0-2.0*nu);
        double g1 = E*(1.0-nu)/TEMP1;
        double g2 = E*nu/TEMP1;
        G[0] = g1; G[7] = g1; G[14] = g1;
        G[1] = g2; G[2] = g2; G[6] = g2; G[8] = g2; G[12] = g2; G[13] = g2;
        G[21] = GG; G[28] = GG; G[35] = GG;
        
        // Build C matrices (4 of 6x3), stored in C(72)
        double C[72];
        for (int i = 0; i < 72; i++) C[i] = 0;
        for (int i = 0; i < 4; i++) {
            int J = 18*i;
            // H is 1-indexed in Fortran: H(I+4), H(I+8), H(I+12)
            // Hinv is 0-indexed here; Fortran H(I+4) with I=1..4 -> H[4..7] 0-indexed
            // Hinv stored col-major: Hinv[row + col*4]
            // Fortran H(I+4) means element at position I+4 in 1-indexed = index I+3 in 0-indexed
            // H is 4x4 col-major: H(1)=Hinv[0], H(2)=Hinv[1], H(3)=Hinv[2], H(4)=Hinv[3]
            // H(5)=Hinv[4], ... H(8)=Hinv[7], H(9)=Hinv[8]...
            // So H(I+4) for I=1..4 -> Hinv[I+3] (0-indexed)
            // H(I+8) -> Hinv[I+7]
            // H(I+12) -> Hinv[I+11]
            double h_i4 = Hinv[i+3];
            double h_i8 = Hinv[i+7];
            double h_i12 = Hinv[i+11];
            // C is 6x3 col-major: C(J+1) means index J (0-indexed)
            // C(J+1)=H(I+4), C(J+5)=H(I+8), C(J+9)=H(I+12)
            // C(J+11)=H(I+12), C(J+12)=H(I+8), C(J+13)=H(I+12)
            // C(J+15)=H(I+4), C(J+16)=H(I+8), C(J+17)=H(I+4)
            C[J+0] = h_i4;
            C[J+4] = h_i8;
            C[J+8] = h_i12;
            C[J+10] = h_i12;
            C[J+11] = h_i8;
            C[J+12] = h_i12;
            C[J+14] = h_i4;
            C[J+15] = h_i8;
            C[J+16] = h_i4;
        }
        
        // HDETER division
        // IOPT is 11..22, so IOPT>=11 && IOPT<=22 -> branch 603
        HDETER = HDETER/36.0;
        if (IOPT <= 16) HDETER = HDETER*2.0;
        
        // KIJ (6x6) zeroed
        double KIJ[36];
        for (int i = 0; i < 36; i++) KIJ[i] = 0;
        
        // Pivot point: NPVT. In test, we need to determine NPVT.
        // The pivot is the first grid point. Let's set NPVT = nodes[0] (grid SIL of node 1)
        // Actually NPVT is set externally. For a standalone test, we need to pick one.
        // The routine inserts stiffness for all 4 points relative to pivot.
        // For the full 18x18 matrix, we need to accumulate all contributions.
        
        // Determine pivot: find which of the 4 nodes is NPVT
        // For our test, let's use NPVT = 1 (grid SIL 1)
        int NPVT = 1;
        
        int NPOINT = -1;
        for (int i = 1; i <= 4; i++) {
            // NECPT(I+1) is grid SIL of node i (1-indexed: NECPT(2..5))
            if (nodes[i-1] == NPVT) {
                NPOINT = 18*i - 35; // 1-indexed position in C
                break;
            }
        }
        if (NPOINT < 0) continue; // pivot not in this tetra
        
        // NPOINT is 1-indexed in C. Convert to 0-indexed: NPOINT-1
        int npoint0 = NPOINT - 1;
        
        // Compute GCT = C(NPOINT) * G (6x3 * 6x6 -> wait)
        // GMMATD(C(NPOINT),6,3,1, G(1),6,6,0, GCT(1))
        // A = C(NPOINT), dims 6x3, transA=1 -> A^T (3x6)
        // B = G, dims 6x6, transB=0 -> G (6x6)
        // result GCT = A^T * G = 3x6
        double GCT[18];
        GMMATD(&C[npoint0], 3, 6, 1, G, 6, 6, 0, GCT);
        
        // GCT *= HDETER
        for (int i = 0; i < 18; i++) GCT[i] *= HDETER;
        
        // Loop through 4 points
        for (int i = 1; i <= 4; i++) {
            // CSID check: NECPT(4*I+3) - for I=1..4, that's NECPT(7,11,15,19)
            // These are CSID-1..4. In our ECPT, CSID for node i is at ECPT[8 + (i-1)*4]
            int csid = (int)ECPT[8 + (i-1)*4];
            double T[9];
            if (csid != 0) {
                TRANSD(csid, T);
                // GMMATD(C(18*I-17),6,3,0, T(1),3,3,0, CT(1))
                // A = C(18*I-17), 6x3, transA=0 -> 6x3
                // B = T, 3x3, transB=0 -> 3x3
                // CT = C * T = 6x3
                double CT[18];
                GMMATD(&C[18*(i-1)], 6, 3, 0, T, 3, 3, 0, CT);
                // GMMATD(GCT(1),3,6,0, CT(1),6,3,0, T(1))
                // A = GCT, 3x6, transA=0 -> 3x6
                // B = CT, 6x3, transB=0 -> 6x3
                // T = GCT * CT = 3x3
                GMMATD(GCT, 3, 6, 0, CT, 6, 3, 0, T);
            } else {
                // GMMATD(GCT(1),3,6,0, C(18*I-17),6,3,0, T(1))
                GMMATD(GCT, 3, 6, 0, &C[18*(i-1)], 6, 3, 0, T);
            }
            
            // Insert 3x3 T into KIJ (6x6)
            // KIJ(1)=T(1), KIJ(2)=T(2), KIJ(3)=T(3)
            // KIJ(7)=T(4), KIJ(8)=T(5), KIJ(9)=T(6)
            // KIJ(13)=T(7), KIJ(14)=T(8), KIJ(15)=T(9)
            // KIJ is 6x6 col-major
            KIJ[0] = T[0]; KIJ[1] = T[1]; KIJ[2] = T[2];
            KIJ[6] = T[3]; KIJ[7] = T[4]; KIJ[8] = T[5];
            KIJ[12] = T[6]; KIJ[13] = T[7]; KIJ[14] = T[8];
            
            // Insert into global K
            // SMA1B(KIJ, NECPT(I+2), -1, IFKGG, 0.0)
            // NECPT(I+2) is grid SIL of node i
            int grid_i = nodes[i-1];
            // Insert 6x6 KIJ into global K at position (grid_i-1)*3
            int base = (grid_i - 1) * 3;
            for (int r = 0; r < 6; r++)
                for (int c = 0; c < 6; c++)
                    K[base+r][base+c] += KIJ[c*6+r]; // KIJ col-major
        }
    }
    
    // Print JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 18; i++) {
        std::cout << "[";
        for (int j = 0; j < 18; j++) {
            std::cout << std::scientific << std::setprecision(6) << K[i][j];
            if (j < 17) std::cout << ",";
        }
        std::cout << "]";
        if (i < 17) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}