#ifndef PML_H
#define PML_H

#include <vector>

// Structure to hold PML parameters
struct PMLParams {
    int thickness;
    double sigma_max;
    double kappa_max;
    double alpha_max; // Use alpha_max for standard CPML
    double alpha_min;
};

class PML {
public:
    // Constructor: sets up all PML coefficients and fields
    PML(int nx, int ny, const PMLParams& params, double dt);

    // Public member variables accessed by update functions
    
    // Scaling factors (1D)
    std::vector<double> kappa_x, kappa_y;

    // CPML update coefficients (1D)
    std::vector<double> b_x, b_y;
    std::vector<double> c_x, c_y;

    // Auxiliary convolution fields (2D)
    std::vector<std::vector<double>> psi_exy, psi_eyx;
    std::vector<std::vector<double>> psi_hzx, psi_hzy;

private:
    // Helper function to generate the PML profiles
    void applyProfile(std::vector<double>& vec, int dimSize, int pmlThickness, double maxVal, double baseVal = 0.0);
};

#endif // PML_H