int main() {
    try {
        // Test case data
        Vector3d nodeA(0.0, 0.0, 0.0);
        Vector3d nodeB(2.0, 0.0, 0.0);
        Real E, G, nu, rho, alpha;
        MAT(1, E, G, nu, rho, alpha); // Get E and G from MAT
        Real A = 0.01;         // m^2
        Real J = 5.0e-6;       // m^4
        
        // Compute stiffness matrix
        Matrix6d K = compute_crod_stiffness(nodeA, nodeB, E, A, G, J);
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 6; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << "[";
            for (int j = 0; j < 6; ++j) {
                if (j > 0) std::cout << ",";
                std::cout << format_double(K(i,j));
            }
            std::cout << "]";
        }
        std::cout << "]}" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}