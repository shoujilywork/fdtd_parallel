#include <iostream>
#include <cmath>
#include <vector>
#include "geometry.h"
#include "materials.h"
#include "sources.h"
#include "utils.h"
#include "pml.h"

// Grid size
const int nx = 800;
const int ny = 400;
const int tmax = 1800;

// Space/time steps
const double dx = 1.0;
const double dt = 0.2 * dx; // 无量纲时间步长

int main() {
    // Initialize field and material grids
    std::vector<std::vector<double>> Ex, Ey, Hz;
    std::vector<std::vector<MaterialType>> materialMap;

    initializeEx(Ex, nx, ny);
    initializeEy(Ey, nx, ny);
    initializeHz(Hz, nx, ny);
    initializeMaterialMap(materialMap, nx, ny);

    // Create two V-shaped antennas
    createVShapedAntenna(nx/2 - 100, ny / 5, 56, 45, 8, true, materialMap);
    createVShapedAntenna(nx/2 + 100, ny / 5, 56, 45, 8, false, materialMap);

    // Save material distribution to VTK file
    saveMaterialMapToVTK(materialMap, 0);

    // Initialize Drude model parameters for copper
    DrudeModel cu_drude;
    cu_drude.omega_p = 1.35e16;
    cu_drude.gamma = 4.1e13;
    cu_drude.eps_inf = 1.0;

    // Compute FDTD update coefficients
    std::vector<std::vector<double>> coeffs = computeUpdateCoefficients(materialMap, cu_drude, dt);

    // Set PML parameters
    PMLParams pmlParams;
    pmlParams.thickness = 4;
    pmlParams.sigma_max = 10.0; // Adjusted for better absorption
    pmlParams.kappa_max = 7.0;
    pmlParams.alpha_min = 0.0; // Alpha is often 0 in many CPML formulations
    
    // Create a complete PML instance
    // NOTE: Your pml.cpp must be updated to initialize all necessary 2D psi fields and coefficients
    PML pml(nx, ny, pmlParams, dt);

    // Initialize source list
    SourceList sources;
    sources.emplace_back(nx/2 - 100, ny / 5+4, 40, 4, 0, 0);
    sources.emplace_back(nx/2 + 100, ny / 5+4, 40, 4, 0, 0);

    // Time iteration loop
    for (int t = 0; t < tmax; ++t) {
        // Pass the entire pml object to the update functions
        updateH(Hz, Ex, Ey, pml, dt, nx, ny);
        updateE(Ex, Ey, Hz, coeffs, pml, materialMap, nx, ny);

        // Apply sources
        if (t < 900) {
            applySources(Ex, Ey, sources, t);
        }

        // Output field snapshots
        if (t % 5 == 0) {
            saveToVTK(Ex,Ey, Hz, t);

        }
    }

    std::cout << "Simulation complete.\n";
    return 0;
}