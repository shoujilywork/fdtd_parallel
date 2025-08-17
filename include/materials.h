#ifndef MATERIALS_H
#define MATERIALS_H

#include <vector>
#include "types.h" 

// 计算 FDTD 更新系数（考虑 Drude 模型）
std::vector<std::vector<double>> computeUpdateCoefficients(
    const std::vector<std::vector<MaterialType>>& materialMap,
    const DrudeModel& drude,
    double dt);

#endif // MATERIALS_H