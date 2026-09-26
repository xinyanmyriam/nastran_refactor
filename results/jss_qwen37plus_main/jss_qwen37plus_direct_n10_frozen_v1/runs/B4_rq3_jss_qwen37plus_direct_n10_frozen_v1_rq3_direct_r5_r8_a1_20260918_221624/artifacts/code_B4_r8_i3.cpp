#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declarations
void ktrmem(int flag);
void sma1b(const Eigen::Matrix<double, 2, 2>& kij, int i_node, int j_node, int npvt, int ifkgg, double temp);

// Global data structures (mimicking Fortran COMMON blocks)
struct CommonData {
    // /CONDAS/
    std::vector<double> consts; // 5 elements

    // /SMA1HT/
    bool heat;

    // /SMA1ET/
    std::vector<double> ecpt; // ECPT(100)

    // /SMA1IO/
    std::vector<double> dum1; // 10 elements
    int ifkgg;
    std::vector<double> dum2; // 1 element
    int if4gg;
    std::vector<double> dum3; // 23 elements

    // /SMA1CL/
    int iopt4;
    int k4ggsw;
    int npvt;
    std::vector<double> dumcl; // 7 elements
    std::vector<int> link; // 10 elements
    int idetck;
    bool dodet;
    int nogO; // renamed to avoid conflict with 'no' keyword

    // /SMA1DP/
    std::vector<double> kij; // KIJ(36)
    std::vector<double> dum7; // 156 elements
    std::vector<double> ksum; // KSUM(36)
    double temp;
    double cosang;
    double sinang;
    double vecl;
    std::vector<double> ivec; // IVEC(3)
    std::vector<double> jvec; // JVEC(3)
    std::vector<double> kvec; // KVEC(3)
    std::vector<double> pvec; // PVEC(3)
    std::vector<double> vsubk; // VSUBK(3)
    std::vector<double> v; // V(3)
    std::vector<double> si; // SI(3)
    int npivot;
    int mpoint;
    int mi;
    int nsubsc;
    std::vector<int> ngrid; // NGRID(4)
    double u1;
    double u2;
    std::vector<double> coord; // COORD(16)
    std::vector<double> dumm8; // 248 elements

    // /MATIN/
    int matid;
    int inflag;
    double eltemp;
    double stress;
    double sinth;
    double costh;

    // /MATOUT/
    std::vector<double> dum99; // 11 elements
    double gsube;
    std::vector<double> dum88; // 6 elements

    // /K3X3/ (equivalenced to kij)
    std::vector<double> k3x3; // 27 elements

    // Constructor to initialize vectors
    CommonData() :
        consts(5, 0.0),
        heat(false),
        ecpt(100, 0.0),
        dum1(10, 0.0),
        ifkgg(0),
        dum2(1, 0.0),
        if4gg(0),
        dum3(23, 0.0),
        iopt4(0),
        k4ggsw(0),
        npvt(0),
        dumcl(7, 0.0),
        link(10, 0),
        idetck(0),
        dodet(false),
        nogO(0),
        kij(36, 0.0),
        dum7(156, 0.0),
        ksum(36, 0.0),
        temp(0.0),
        cosang(0.0),
        sinang(0.0),
        vecl(0.0),
        ivec(3, 0.0),
        jvec(3, 0.0),
        kvec(3, 0.0),
        pvec(3, 0.0),
        vsubk(3, 0.0),
        v(3, 0.0),
        si(3, 0.0),
        npivot(0),
        mpoint(0),
        mi(0),
        nsubsc(0),
        ngrid(4, 0),
        u1(0.0),
        u2(0.0),
        coord(16, 0.0),
        dumm8(248, 0.0),
        matid(0),
        inflag(0),
        eltemp(0.0),
        stress(0.0),
        sinth(0.0),
        costh(0.0),
        dum99(11, 0.0),
        gsube(0.0),
        dum88(6, 0.0),
        k3x3(27, 0.0) {}
};

CommonData common;

// Helper function to get DEGRA (degrees to radians conversion factor)
double get_degra() {
    return M_PI / 180.0;
}

// Helper function to set up the mapping matrix M (12 elements)
std::vector<int> get_m() {
    return {1, 2, 4, 2, 3, 1, 3, 4, 2, 4, 1, 3};
}

// Compute triangle normal and orthonormal basis (s,t,n)
std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
compute_triangle_basis(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
    Eigen::Vector3d e1 = p1 - p0;
    Eigen::Vector3d e2 = p2 - p0;
    Eigen::Vector3d n = e1.cross(e2);
    double norm_n = n.norm();
    if (norm_n < 1e-12) {
        // Degenerate triangle
        return {Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(), Eigen::Vector3d::UnitZ()};
    }
    n /= norm_n;

    Eigen::Vector3d s = e1.normalized();
    Eigen::Vector3d t = n.cross(s);
    t.normalize();

    return {s, t, n};
}

// Compute 6x6 CST stiffness in global coordinates, then extract 2x2 blocks
// Returns vector of nine 2x2 blocks: [k11, k12, k13, k21, k22, k23, k31, k32, k33]
std::vector<Eigen::Matrix<double, 2, 2>> compute_triangular_stiffness_blocks(
    const std::vector<Eigen::Vector3d>& coords,
    double e, double nu, double t) {

    // coords has 3 points: p0, p1, p2
    const auto& p0 = coords[0];
    const auto& p1 = coords[1];
    const auto& p2 = coords[2];

    // Area
    double area = 0.5 * (p1 - p0).cross(p2 - p0).norm();
    if (area < 1e-15) {
        std::vector<Eigen::Matrix<double, 2, 2>> blocks(9, Eigen::Matrix<double, 2, 2>::Zero());
        return blocks;
    }

    // Material matrix D (plane stress)
    double d11 = e / (1.0 - nu*nu);
    double d12 = d11 * nu;
    double d22 = d11;
    double d66 = e / (2.0 * (1.0 + nu));

    Eigen::Matrix<double, 3, 3> d;
    d << d11, d12, 0.0,
         d12, d22, 0.0,
         0.0, 0.0, d66;

    // Strain-displacement B matrix (3x6)
    double x1 = p0(0), y1 = p0(1);
    double x2 = p1(0), y2 = p1(1);
    double x3 = p2(0), y3 = p2(1);

    double b1 = y2 - y3;
    double b2 = y3 - y1;
    double b3 = y1 - y2;
    double c1 = x3 - x2;
    double c2 = x1 - x3;
    double c3 = x2 - x1;

    Eigen::Matrix<double, 3, 6> b;
    b << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;
    b /= (2.0 * area);

    // Full 6x6 stiffness
    Eigen::Matrix<double, 6, 6> k_full = t * area * b.transpose() * d * b;

    // Extract 2x2 blocks: kIJ for node I,J ∈ {0,1,2}
    std::vector<Eigen::Matrix<double, 2, 2>> blocks;
    for (int I = 0; I < 3; ++I) {
        for (int J = 0; J < 3; ++J) {
            Eigen::Matrix<double, 2, 2> blk;
            blk << k_full(2*I+0, 2*J+0), k_full(2*I+0, 2*J+1),
                   k_full(2*I+1, 2*J+0), k_full(2*I+1, 2*J+1);
            blocks.push_back(blk);
        }
    }
    return blocks;
}

// Compute 3x3 rotation matrix T: [s_x s_y 0; t_x t_y 0; 0 0 1]
Eigen::Matrix<double, 3, 3> compute_T_matrix(const Eigen::Vector3d& s, const Eigen::Vector3d& t) {
    Eigen::Matrix<double, 3, 3> T;
    T << s(0), s(1), 0.0,
         t(0), t(1), 0.0,
         0.0,  0.0,  1.0;
    return T;
}

// Rotate a 2x2 block: k_global = T2x2^T * k_local * T2x2, where T2x2 = [[sx,sy],[tx,ty]]
Eigen::Matrix<double, 2, 2> rotate_2x2_block(const Eigen::Matrix<double, 2, 2>& k_local,
                                              double sx, double sy, double tx, double ty) {
    Eigen::Matrix<double, 2, 2> T2x2;
    T2x2 << sx, sy,
            tx, ty;
    return T2x2.transpose() * k_local * T2x2;
}

// Mock implementation of KTRMEM - computes triangular membrane stiffness
void ktrmem(int flag) {
    // Get the 3 node coordinates from ECPT
    std::vector<Eigen::Vector3d> coords;
    for (int i = 0; i < 3; ++i) {
        double x = common.ecpt[4*i + 10]; // ECPT(11), ECPT(15), ECPT(19)
        double y = common.ecpt[4*i + 11]; // ECPT(12), ECPT(16), ECPT(20)
        double z = common.ecpt[4*i + 12]; // ECPT(13), ECPT(17), ECPT(21)
        coords.emplace_back(x, y, z);
    }

    // Get material properties
    double e = common.consts[0]; // E
    double nu = common.consts[1]; // nu
    double t = common.ecpt[7]; // thickness

    // Compute local stiffness blocks (6x6 → nine 2x2)
    auto blocks = compute_triangular_stiffness_blocks(coords, e, nu, t);

    // Compute triangle basis: s, t, n
    auto [s, t_vec, n] = compute_triangle_basis(coords[0], coords[1], coords[2]);
    Eigen::Matrix<double, 3, 3> T = compute_T_matrix(s, t_vec);
    double sx = s(0), sy = s(1), tx = t_vec(0), ty = t_vec(1);

    // Rotate each 2x2 block and store in k3x3 (27 elements: 9 blocks × 2×2 = 36? but Fortran uses 27 → we store only upper 3×3 of 6x6? No — per spec: 27 = 9 blocks × 3 entries? Actually: legacy stores 3×3 per block, but we need 2×2. So reinterpret k3x3 as 9 × (2×2) = 36 → but size is 27. So Fortran packs as 3×3 per block, but only uses first 2 rows/cols. We'll store 2x2 in first 4 entries of each 3x3 block, zero-padded.)
    // Instead: use k3x3 as flat array of 27 doubles → store each 2x2 block row-wise in first 4 slots of each 3-element group? No — 27 / 9 = 3 → so each block is 3 entries → implies they store only diagonal? That contradicts spec.

    // Per problem statement: "K3X3(27) which is 3 matrices of 3x3 each" → actually 3 blocks, not 9. But CQUAD4 uses 4 triangles, each contributing 3 blocks? Let's re-read: "The 3x3 T-matrix rotation T^T * K_block * T must be applied to EACH 3x3 block" — and "K3X3(27)" → 27 = 3 × 3 × 3 → three 3×3 blocks.

    // So: each triangle contributes exactly THREE 3×3 blocks: likely k11, k12, k13 (pivot-to-all). That matches pivot-based decomposition.

    // Therefore: for triangle with nodes [P, A, B] (P=pivot), compute kPP, kPA, kPB — three 3×3 blocks.

    // But our CST gives 6×6 → k11, k12, k13, k21, k22, k23, k31, k32, k33. For pivot-based, we need k11, k12, k13 (if node 0 is pivot).

    // So assume coords[0] = pivot, coords[1] = node A, coords[2] = node B.

    // Extract k11, k12, k13 (each 2x2) from k_full
    std::vector<Eigen::Matrix<double, 2, 2>> target_blocks = {
        blocks[0], // k11
        blocks[1], // k12
        blocks[2]  // k13
    };

    // Rotate each 2x2 block
    std::vector<Eigen::Matrix<double, 2, 2>> rotated_blocks;
    for (const auto& blk : target_blocks) {
        rotated_blocks.push_back(rotate_2x2_block(blk, sx, sy, tx, ty));
    }

    // Store in k3x3: each 2x2 block stored as 4 values → but k3x3 has 27 entries → we only need 3 blocks × 4 = 12. So pad with zeros or use first 12.
    // However, Fortran expects 27, and legacy code indexes k3x3(1..27). So store each 2x2 as:
    // k3x3[0..3] = k11 (row0: k11_00,k11_01,0; row1: k11_10,k11_11,0; row2: 0,0,0) → no, too wasteful.

    // Instead: reinterpret — the problem says "3x3 T-matrix rotation", so they *do* use 3×3 blocks. But membrane has only 2 DOF/node → third DOF is unused (drill? no — CQUAD4 membrane has 2 DOF). So third row/col is zero.

    // So build 3×3 block from 2×2:
    auto to_3x3 = [](const Eigen::Matrix<double, 2, 2>& blk) -> Eigen::Matrix<double, 3, 3> {
        Eigen::Matrix<double, 3, 3> m = Eigen::Matrix<double, 3, 3>::Zero();
        m.block<2,2>(0,0) = blk;
        return m;
    };

    // Now rotate full 3x3 blocks
    for (int idx = 0; idx < 3; ++idx) {
        Eigen::Matrix<double, 3, 3> k33_local = to_3x3(rotated_blocks[idx]);
        // Apply T^T * k33_local * T
        Eigen::Matrix<double, 3, 3> k33_global = T.transpose() * k33_local * T;

        // Store in k3x3: block idx starts at idx*9
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                common.k3x3[idx*9 + i*3 + j] = k33_global(i, j);
            }
        }
    }
}

// Accumulate stiffness into global 8x8 matrix
void sma1b(const Eigen::Matrix<double, 2, 2>& kij, int i_node, int j_node, int npvt, int ifkgg, double temp) {
    static Eigen::Matrix<double, 8, 8> global_k = Eigen::Matrix<double, 8, 8>::Zero();

    // Map node ID (1-based) to DOF indices: node i → DOFs (2*i-2, 2*i-1) 0-based
    int dofi0 = (i_node - 1) * 2;
    int dofi1 = dofi0 + 1;
    int dofj0 = (j_node - 1) * 2;
    int dofj1 = dofj0 + 1;

    // Add 2x2 block
    global_k(dofi0, dofj0) += kij(0,0);
    global_k(dofi0, dofj1) += kij(0,1);
    global_k(dofi1, dofj0) += kij(1,0);
    global_k(dofi1, dofj1) += kij(1,1);

    // Also add transpose if i != j (symmetric)
    if (i_node != j_node) {
        global_k(dofj0, dofi0) += kij(0,0);
        global_k(dofj0, dofi1) += kij(1,0);
        global_k(dofj1, dofi0) += kij(0,1);
        global_k(dofj1, dofi1) += kij(1,1);
    }
}

// Main CQUAD4 stiffness computation function
Eigen::Matrix<double, 8, 8> compute_cquad4_stiffness(
    const std::vector<std::vector<double>>& nodes,
    double e, double nu, double t) {

    // Initialize common data
    common.consts[0] = e;
    common.consts[1] = nu;
    common.consts[3] = get_degra(); // DEGRA

    // Set up ECPT array with test data
    common.ecpt[0] = 1.0; // EL.ID
    common.ecpt[1] = 1.0; // A
    common.ecpt[2] = 2.0; // B
    common.ecpt[3] = 3.0; // C
    common.ecpt[4] = 4.0; // D
    common.ecpt[5] = 0.0; // THETA
    common.ecpt[6] = 1.0; // MATID
    common.ecpt[7] = t;   // T
    common.ecpt[8] = 0.0; // NON-STRUCT. MASS
    common.ecpt[9] = 0.0; // COORD SYS ID

    // Coordinates: nodes are (0,0), (2,0), (2,1.5), (0,1.5)
    for (int i = 0; i < 4; ++i) {
        common.ecpt[10 + 4*i + 0] = nodes[i][0]; // X
        common.ecpt[10 + 4*i + 1] = nodes[i][1]; // Y
        common.ecpt[10 + 4*i + 2] = 0.0;         // Z
    }

    // Set other required values
    common.npvt = 1; // pivot point = node 1
    common.iopt4 = 0;
    common.gsube = 0.0;

    // Copy coordinates to COORD array
    for (int i = 0; i < 4; ++i) {
        common.coord[4*i + 0] = static_cast<double>(i+1); // CSID
        common.coord[4*i + 1] = nodes[i][0]; // X
        common.coord[4*i + 2] = nodes[i][1]; // Y
        common.coord[4*i + 3] = 0.0;         // Z
    }

    // Set up NGRID
    common.ngrid[0] = 1;
    common.ngrid[1] = 2;
    common.ngrid[2] = 3;
    common.ngrid[3] = 4;

    // HRING = false
    bool hring = false;

    // Compute angle and direction vectors
    double angl = common.ecpt[5] * common.consts[3]; // THETA * DEGRA
    common.cosang = std::cos(angl);
    common.sinang = std::sin(angl);

    // IVEC = ECPT(15)-ECPT(11), ECPT(16)-ECPT(12), ECPT(17)-ECPT(13)
    common.ivec[0] = common.ecpt[13] - common.ecpt[10]; // X2-X1
    common.ivec[1] = common.ecpt[14] - common.ecpt[11]; // Y2-Y1
    common.ivec[2] = common.ecpt[15] - common.ecpt[12]; // Z2-Z1

    common.vecl = std::sqrt(common.ivec[0]*common.ivec[0] +
                           common.ivec[1]*common.ivec[1] +
                           common.ivec[2]*common.ivec[2]);
    if (common.vecl != 0.0) {
        common.ivec[0] /= common.vecl;
        common.ivec[1] /= common.vecl;
        common.ivec[2] /= common.vecl;
    }

    // VSUBK = IVEC × (ECPT(25)-ECPT(13), ECPT(24)-ECPT(12), ECPT(23)-ECPT(11))
    // ECPT(23-25) = X4,Y4,Z4
    common.vsubk[0] = common.ivec[1]*(common.ecpt[24]-common.ecpt[12]) -
                      common.ivec[2]*(common.ecpt[23]-common.ecpt[11]);
    common.vsubk[1] = common.ivec[2]*(common.ecpt[22]-common.ecpt[10]) -
                      common.ivec[0]*(common.ecpt[24]-common.ecpt[12]);
    common.vsubk[2] = common.ivec[0]*(common.ecpt[23]-common.ecpt[11]) -
                      common.ivec[1]*(common.ecpt[22]-common.ecpt[10]);

    common.vecl = std::sqrt(common.vsubk[0]*common.vsubk[0] +
                           common.vsubk[1]*common.vsubk[1] +
                           common.vsubk[2]*common.vsubk[2]);
    if (common.vecl != 0.0) {
        common.kvec[0] = common.vsubk[0] / common.vecl;
        common.kvec[1] = common.vsubk[1] / common.vecl;
        common.kvec[2] = common.vsubk[2] / common.vecl;
    }

    // JVEC = KVEC × IVEC
    common.jvec[0] = common.kvec[1]*common.ivec[2] - common.kvec[2]*common.ivec[1];
    common.jvec[1] = common.kvec[2]*common.ivec[0] - common.kvec[0]*common.ivec[2];
    common.jvec[2] = common.kvec[0]*common.ivec[1] - common.kvec[1]*common.ivec[0];

    // PVEC = COSANG*IVEC + SINANG*JVEC
    for (int i = 0; i < 3; ++i) {
        common.pvec[i] = common.cosang * common.ivec[i] + common.sinang * common.jvec[i];
    }

    // Set up ECPT for membrane use
    common.ecpt[4] = common.ecpt[5];
    common.ecpt[5] = common.ecpt[6];
    if (!hring) {
        common.ecpt[6] = common.ecpt[7] / 2.0;
    }
    common.ecpt[7] = common.ecpt[8];

    // Find pivot point
    int npivot = 0;
    for (int i = 0; i < 4; ++i) {
        if (common.npvt == common.ngrid[i]) {
            npivot = i + 1; // 1-based
            break;
        }
    }

    // Compute JNOT
    int jnot = 0;
    if (npivot <= 2) {
        jnot = npivot + 2;
    } else {
        jnot = npivot - 2;
    }

    // Zero out KSUM
    for (int i = 0; i < 36; ++i) {
        common.ksum[i] = 0.0;
    }

    // Get mapping matrix
    std::vector<int> m = get_m();

    // Reset global stiffness
    static Eigen::Matrix<double, 8, 8> global_k = Eigen::Matrix<double, 8, 8>::Zero();
    global_k.setZero();

    // Loop through 4 triangles: pivot + two others
    // Triangles: (1,2,3), (1,3,4), (2,3,4), (1,2,4) — but standard is pivot + consecutive pairs
    // Legacy uses: for pivot=1: triangles (1,2,3), (1,3,4), (1,4,2) — but jnot skips one.
    // Instead, use canonical four: (1,2,3), (1,3,4), (2,3,4), (1,2,4) — but skip jnot.

    // Build triangle node lists
    std::vector<std::vector<int>> triangles = {
        {1,2,3}, // tri1
        {1,3,4}, // tri2
        {2,3,4}, // tri3
        {1,2,4}  // tri4
    };

    for (int tri_idx = 0; tri_idx < 4; ++tri_idx) {
        if (tri_idx + 1 == jnot) continue;

        const auto& tri = triangles[tri_idx];
        // Fill ECPT for triangle: coords of tri[0], tri[1], tri[2]
        for (int i = 0; i < 3; ++i) {
            int node_id = tri[i];
            // node_id is 1-based; get coords from nodes[node_id-1]
            int idx_in_nodes = node_id - 1;
            common.ecpt[10 + 4*i + 0] = nodes[idx_in_nodes][0];
            common.ecpt[10 + 4*i + 1] = nodes[idx_in_nodes][1];
            common.ecpt[10 + 4*i + 2] = 0.0;
        }

        // Call KTRMEM
        ktrmem(1);

        // Extract and assemble three 2x2 blocks: k11, k12, k13 (pivot is first node: tri[0])
        int pivot_node = tri[0];
        int nodeA = tri[1];
        int nodeB = tri[2];

        // Block 0: k11 (pivot-pivot) → 2x2
        Eigen::Matrix<double, 2, 2> k11;
        k11 << common.k3x3[0], common.k3x3[1],
               common.k3x3[3], common.k3x3[4];
        sma1b(k11, pivot_node, pivot_node, common.npvt, common.ifkgg, common.temp);

        // Block 1: k12 (pivot-nodeA) → 2x2
        Eigen::Matrix<double, 2, 2> k12;
        k12 << common.k3x3[9], common.k3x3[10],
               common.k3x3[12], common.k3x3[13];
        sma1b(k12, pivot_node, nodeA, common.npvt, common.ifkgg, common.temp);

        // Block 2: k13 (pivot-nodeB) → 2x2
        Eigen::Matrix<double, 2, 2> k13;
        k13 << common.k3x3[18], common.k3x3[19],
               common.k3x3[21], common.k3x3[22];
        sma1b(k13, pivot_node, nodeB, common.npvt, common.ifkgg, common.temp);
    }

    return global_k;
}

int main() {
    // Test case: Nodes (0,0), (2,0), (2,1.5), (0,1.5)
    std::vector<std::vector<double>> nodes = {
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.5},
        {0.0, 1.5}
    };

    double e = 200e9;      // Pa
    double nu = 0.3;       // Poisson's ratio
    double t = 0.01;       // m

    // Compute stiffness matrix
    Eigen::Matrix<double, 8, 8> k = compute_cquad4_stiffness(nodes, e, nu, t);

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 8; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 8; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(15) << k(i,j);
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}