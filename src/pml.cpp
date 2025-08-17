#include "pml.h"
#include <cmath>
#include <vector>
#include <iostream>

/**
 * @brief PML Constructor.
 * This is where all the magic happens. It calculates the PML profiles for sigma, kappa,
 * and alpha, and then uses them to pre-calculate the update coefficients 'b' and 'c'
 * for the CPML formulation. It also correctly resizes the 2D psi auxiliary fields.
 */
PML::PML(int nx, int ny, const PMLParams& params, double dt) {
    // --- 1. Resize all member vectors to their correct dimensions ---
    kappa_x.resize(nx); kappa_y.resize(ny);
    b_x.resize(nx);     b_y.resize(ny);
    c_x.resize(nx);     c_y.resize(ny);

    // Resize 2D psi fields and initialize all values to 0.0
    psi_exy.assign(nx, std::vector<double>(ny, 0.0));
    psi_eyx.assign(nx, std::vector<double>(ny, 0.0));
    psi_hzx.assign(nx, std::vector<double>(ny, 0.0));
    psi_hzy.assign(nx, std::vector<double>(ny, 0.0));

    // --- 2. Create temporary profile vectors for sigma and alpha ---
    std::vector<double> sigma_x(nx), sigma_y(ny);
    std::vector<double> alpha_x(nx), alpha_y(ny);

    // --- 3. Generate the PML profiles using the helper function ---
    // Base value is 0.0 for sigma and alpha, 1.0 for kappa
    applyProfile(sigma_x, nx, params.thickness, params.sigma_max, 0.0);
    applyProfile(sigma_y, ny, params.thickness, params.sigma_max, 0.0);
    applyProfile(kappa_x, nx, params.thickness, params.kappa_max, 1.0);
    applyProfile(kappa_y, ny, params.thickness, params.kappa_max, 1.0);
    applyProfile(alpha_x, nx, params.thickness, params.alpha_max, 0.0);
    applyProfile(alpha_y, ny, params.thickness, params.alpha_max, 0.0);

    // --- 4. Calculate the 'b' and 'c' CPML coefficients ---
    // For the x-direction
    for (int i = 0; i < nx; ++i) {
        // In normalized units, we assume epsilon_0 = 1
        double sigma_over_kappa = sigma_x[i] / kappa_x[i];
        b_x[i] = exp(-(sigma_over_kappa + alpha_x[i]) * dt);
        
        // Denominator for the 'c' coefficient calculation
        double denom = sigma_x[i] * kappa_x[i] + alpha_x[i] * kappa_x[i] * kappa_x[i];
        
        // Avoid division by zero in the main grid where sigma and alpha are 0
        if (std::abs(denom) < 1e-9) {
            c_x[i] = 0.0;
        } else {
            c_x[i] = sigma_x[i] * (b_x[i] - 1.0) / denom;
        }
    }

    // For the y-direction
    for (int j = 0; j < ny; ++j) {
        double sigma_over_kappa = sigma_y[j] / kappa_y[j];
        b_y[j] = exp(-(sigma_over_kappa + alpha_y[j]) * dt);

        double denom = sigma_y[j] * kappa_y[j] + alpha_y[j] * kappa_y[j] * kappa_y[j];

        if (std::abs(denom) < 1e-9) {
            c_y[j] = 0.0;
        } else {
            c_y[j] = sigma_y[j] * (b_y[j] - 1.0) / denom;
        }
    }
    std::cout << "PML initialized successfully." << std::endl;
}

/**
 * @brief Generic helper function to create a cubic PML parameter profile.
 * The value is 'baseVal' in the main grid and ramps up cubically inside the PML regions.
 * @param vec The vector to fill (e.g., sigma_x, kappa_y).
 * @param dimSize The size of the dimension (nx or ny).
 * @param pmlThickness The thickness of the PML in grid cells.
 * @param maxVal The maximum value of the parameter at the simulation edge.
 * @param baseVal The value in the non-PML region (e.g., 0.0 for sigma, 1.0 for kappa).
 */
void PML::applyProfile(std::vector<double>& vec, int dimSize, int pmlThickness, double maxVal, double baseVal) {
    for (int i = 0; i < dimSize; ++i) {
        double pos_in_pml = 0.0;
        
        if (i < pmlThickness) { // Left/Bottom PML region
            pos_in_pml = static_cast<double>(pmlThickness - i) / pmlThickness;
        } else if (i >= dimSize - pmlThickness) { // Right/Top PML region
            pos_in_pml = static_cast<double>(i - (dimSize - pmlThickness - 1)) / pmlThickness;
        }

        if (pos_in_pml > 0.0) {
            vec[i] = baseVal + (maxVal - baseVal) * pow(pos_in_pml, 3);
        } else {
            vec[i] = baseVal;
        }
    }
}