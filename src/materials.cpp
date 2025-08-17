#include "materials.h"
#include <cmath>
#include "utils.h"

std::vector<std::vector<double>> computeUpdateCoefficients(
    const std::vector<std::vector<MaterialType>>& materialMap,
    const DrudeModel& drude,
    double dt) {

    int nx = materialMap.size();
    int ny = materialMap[0].size();
    std::vector<std::vector<double>> coeffs(nx, std::vector<double>(ny));

    // ... (Drude model parameter setup is unchanged) ...

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (materialMap[i][j] == MaterialType::METAL) {
                // The PEC boundary condition is handled in updateE by forcing E=0.
                // Setting the coefficient to 0 here is also a valid way to do it.
                coeffs[i][j] = 0.0;
            } else {
                // --- THIS IS THE FIX ---
                // The update coefficient for vacuum (where permittivity is 1) is simply the time step dt.
                coeffs[i][j] = dt;
            }
        }
    }
    return coeffs;
}