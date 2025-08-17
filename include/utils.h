#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include "types.h"
#include "pml.h"
// 更新 H 场
void updateH(std::vector<std::vector<double>>& Hz,
             const std::vector<std::vector<double>>& Ex,
             const std::vector<std::vector<double>>& Ey,
             PML& pml,
             double dt, // Added dt as an argument
             int nx, int ny);

// 更新 E 场
void updateE(std::vector<std::vector<double>>& Ex,
             std::vector<std::vector<double>>& Ey,
             const std::vector<std::vector<double>>& Hz,
             const std::vector<std::vector<double>>& coeffs,
             PML& pml,
             const std::vector<std::vector<MaterialType>>& materialMap, // Added materialMap
             int nx, int ny);

// 计算 FDTD 更新系数（基于 Drude 模型）
std::vector<std::vector<double>> computeUpdateCoefficients(
    const std::vector<std::vector<MaterialType>>& materialMap,
    const DrudeModel& drude,
    double dt);

// 多激励源支持
void applySources(std::vector<std::vector<double>>& Ex,
                  const std::vector<struct GaussianSource>& sources,
                  int t);

void saveMaterialMapToVTK(const std::vector<std::vector<MaterialType>>& materialMap, int step);

// 输出电场分布到 VTK 文件
void saveToVTK(const std::vector<std::vector<double>>& Ex,
                     const std::vector<std::vector<double>>& Ey,
                     const std::vector<std::vector<double>>& Hz,
                     int step);

// 初始化场
void initializeEx(std::vector<std::vector<double>>& Ex, int nx, int ny);
void initializeEy(std::vector<std::vector<double>>& Ey, int nx, int ny);
void initializeHz(std::vector<std::vector<double>>& Hz, int nx, int ny);
void initializeMaterialMap(std::vector<std::vector<MaterialType>>& materialMap, int nx, int ny);

#endif // UTILS_H