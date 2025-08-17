#include "sources.h"
#include <cmath>
#include <vector>

// Gaussian pulse function (this remains unchanged)
double GaussianSource::pulse(int t) const {
    double tau = t - t0;
    double envelope = exp(-pow(tau, 2) / (2.0 * spread * spread));
    return envelope * cos(phase);
}

/**
 * @brief Applies all sources to both Ex and Ey fields based on their polarization.
 * @param Ex The Ex field grid.
 * @param Ey The Ey field grid.
 * @param sources The list of sources to apply.
 * @param t The current time step.
 */
void applySources(std::vector<std::vector<double>>& Ex,
                  std::vector<std::vector<double>>& Ey,
                  const SourceList& sources,
                  int t) {
    for (const auto& source : sources) {
        if (source.x >= 0 && source.x < (int)Ex.size() &&
            source.y >= 0 && source.y < (int)Ex[0].size()) {
            
            // Calculate the total pulse amplitude at this time step
            double pulse_amplitude = source.pulse(t);

            // Decompose the amplitude into Ex and Ey components using the polarization angle
            double ex_component = pulse_amplitude * cos(source.polarization_angle);
            double ey_component = pulse_amplitude * sin(source.polarization_angle);

            // Add the components to their respective fields
            Ex[source.x][source.y] += ex_component;
            Ey[source.x][source.y] += ey_component;
        }
    }
}