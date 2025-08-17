#ifndef SOURCES_H
#define SOURCES_H

#include <vector>

class GaussianSource {
public:
    int x, y;                // Position on the grid
    double t0;               // Time delay (center of the pulse)
    double spread;           // Width of the pulse
    double phase;            // Phase of the pulse
    double polarization_angle; // Angle of E-field in radians (0=Ex, PI/2=Ey)

    // Constructor
    GaussianSource(int x_pos, int y_pos, double delay, double pulse_width, double ph, double pol_angle)
        : x(x_pos), y(y_pos), t0(delay), spread(pulse_width), phase(ph), polarization_angle(pol_angle) {}

    // Method to get the pulse value at a given time
    double pulse(int t) const;
};

// A list of sources
using SourceList = std::vector<GaussianSource>;

// Function to apply all sources to the E-field grids
void applySources(std::vector<std::vector<double>>& Ex,
                  std::vector<std::vector<double>>& Ey,
                  const SourceList& sources,
                  int t);

#endif // SOURCES_H