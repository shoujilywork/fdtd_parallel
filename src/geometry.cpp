#include "geometry.h"
#include <cmath>
#include "types.h"
extern std::vector<std::vector<MaterialType>> materialMap;
const int nx = 800;  // 网格大小（应与 main.cu 一致）
const int ny = 400;

void createVShapedAntenna(int centerX, int centerY, int armLength,
                          int angleDeg, int width, bool leftFacing,
                          std::vector<std::vector<MaterialType>>& materialMap) {
    double angleRad = angleDeg * M_PI / 180.0;
    double cosA = cos(angleRad / 2);
    double sinA = sin(angleRad / 2);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double dx = i - centerX;
            double dy = j - centerY;

            float proj1 = dx * cosA + dy * sinA;
            float perp1 = -dx * sinA + dy * cosA;

            float proj2 = dx * cosA - dy * sinA;
            float perp2 = dx * sinA + dy * cosA;

            if (leftFacing) {
                if (fabs(proj2) <= armLength && fabs(perp2) <= width / 2) {
                    materialMap[i][j] = MaterialType::METAL;
                }
            } else {
                if (fabs(proj1) <= armLength && fabs(perp1) <= width / 2) {
                    materialMap[i][j] = MaterialType::METAL;
                }
            }
        }
    }
}