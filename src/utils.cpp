#include "utils.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "types.h"
#include "pml.h" // Include PML header to use the class

// Update H field (Hz) with full 2D CPML
void updateH(std::vector<std::vector<double>>& Hz,
             const std::vector<std::vector<double>>& Ex,
             const std::vector<std::vector<double>>& Ey,
             PML& pml,
             double dt, // Added dt as an argument
             int nx, int ny) {
    // 并行化外层循环 (i 循环)
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < nx - 1; ++i) {
        for (int j = 0; j < ny - 1; ++j) {
            double dEy_dx = Ey[i + 1][j] - Ey[i][j];
            double dEx_dy = Ex[i][j + 1] - Ex[i][j];

            pml.psi_hzx[i][j] = pml.b_x[i] * pml.psi_hzx[i][j] + pml.c_x[i] * dEy_dx;
            pml.psi_hzy[i][j] = pml.b_y[j] * pml.psi_hzy[i][j] + pml.c_y[j] * dEx_dy;
            
            // --- FIX: The entire update is now correctly scaled by dt ---
            double curl_E_pml = ( (dEy_dx / pml.kappa_x[i] + pml.psi_hzx[i][j]) - 
                                  (dEx_dy / pml.kappa_y[j] + pml.psi_hzy[i][j]) );
            
            Hz[i][j] -= dt * curl_E_pml; // dt scaling factor applied here
        }
    }
}

/**
 * @brief Combined Update for E fields (Ex, Ey) with full 2D CPML and PEC boundary conditions.
 * @param materialMap The map defining the location of metal and other materials.
 */
void updateE(std::vector<std::vector<double>>& Ex,
             std::vector<std::vector<double>>& Ey,
             const std::vector<std::vector<double>>& Hz,
             const std::vector<std::vector<double>>& coeffs,
             PML& pml,
             const std::vector<std::vector<MaterialType>>& materialMap, // Added materialMap
             int nx, int ny) {

    // --- Step 1: Update E-fields everywhere using the standard update equations ---
    // Update Ex
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            double dHz_dy = Hz[i][j] - Hz[i][j - 1];
            pml.psi_exy[i][j] = pml.b_y[j] * pml.psi_exy[i][j] + pml.c_y[j] * dHz_dy;
            Ex[i][j] += coeffs[i][j] * (dHz_dy / pml.kappa_y[j] + pml.psi_exy[i][j]);
        }
    }

    // Update Ey
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            double dHz_dx = Hz[i][j] - Hz[i - 1][j];
            pml.psi_eyx[i][j] = pml.b_x[i] * pml.psi_eyx[i][j] + pml.c_x[i] * dHz_dx;
            Ey[i][j] -= coeffs[i][j] * (dHz_dx / pml.kappa_x[i] + pml.psi_eyx[i][j]);
        }
    }

    // --- Step 2: 并行施加 PEC 边界条件 ---
    // 由于这是对金属区域的简单赋值，也可以并行化
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (materialMap[i][j] == MaterialType::METAL) {
                Ex[i][j] = 0.0;
                Ey[i][j] = 0.0;
            }
        }
    }
}


// Initialization and VTK save functions remain the same...

// Initialize Ex field to zero
void initializeEx(std::vector<std::vector<double>>& Ex, int nx, int ny) {
    Ex.assign(nx, std::vector<double>(ny, 0.0));
    std::cout << "Ex initialized with size " << nx << "x" << ny << std::endl;
}

// Initialize Ey field to zero
void initializeEy(std::vector<std::vector<double>>& Ey, int nx, int ny) {
    Ey.assign(nx, std::vector<double>(ny, 0.0));
    std::cout << "Ey initialized with size " << nx << "x" << ny << std::endl;
}

// Initialize Hz field to zero
void initializeHz(std::vector<std::vector<double>>& Hz, int nx, int ny) {
    Hz.assign(nx, std::vector<double>(ny, 0.0));
    std::cout << "Hz initialized with size " << nx << "x" << ny << std::endl;
}

// Initialize material grid
void initializeMaterialMap(std::vector<std::vector<MaterialType>>& materialMap, int nx, int ny) {
    materialMap.assign(nx, std::vector<MaterialType>(ny, MaterialType::DIELECTRIC));
    std::cout << "Material map initialized with size " << nx << "x" << ny << std::endl;
}

/**
 * @brief Saves the electric field (Ex, Ey) as vector data to a single VTK file.
 * This allows visualization of the E-field direction and magnitude in tools like ParaView.
 * @param Ex The x-component of the electric field.
 * @param Ey The y-component of the electric field.
 * @param step The current simulation time step, used for the filename.
 */
void saveToVTK(const std::vector<std::vector<double>>& Ex,
                     const std::vector<std::vector<double>>& Ey,
                     const std::vector<std::vector<double>>& Hz,
                     int step) {
    // Create a filename for the current step
    std::string filename = "fdtd_fields_" + std::to_string(step) + ".vtk";
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error: Failed to open output file " << filename << "!" << std::endl;
        return;
    }

    int nx = Ex.size();
    int ny = Ex[0].size();

    // --- VTK Header ---
    fout << "# vtk DataFile Version 3.0\n";
    fout << "2D EM Field Data (E-Vector, H-Scalar)\n";
    fout << "ASCII\n";
    fout << "DATASET STRUCTURED_GRID\n";
    fout << "DIMENSIONS " << nx << " " << ny << " 1\n";
    
    // --- Grid Point Coordinates ---
    fout << "POINTS " << nx * ny << " float\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fout << static_cast<float>(i) << " " << static_cast<float>(j) << " 0.0\n";
        }
    }

    // --- Field Data Section ---
    fout << "POINT_DATA " << nx * ny << "\n";
    
    // --- 1. Vector Field Data for E-field ---
    fout << "VECTORS E_field float\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fout << static_cast<float>(Ex[i][j]) << " " 
                 << static_cast<float>(Ey[i][j]) << " 0.0\n";
        }
    }

    // --- 2. Scalar Field Data for H-field ---
    fout << "SCALARS Hz_field float 1\n"; // "1" specifies one component (scalar)
    fout << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fout << static_cast<float>(Hz[i][j]) << "\n";
        }
    }

    fout.close();
    std::cout << "Saved E/H field data to " << filename << std::endl;
}


// 将 materialMap 保存为 VTK 文件
void saveMaterialMapToVTK(const std::vector<std::vector<MaterialType>>& materialMap, int step) {
    std::string filename = "material_" + std::to_string(step) + ".vtk";
    std::ofstream fout(filename);

    if (!fout) {
        std::cerr << "Error: Cannot open file " << filename << " for writing." << std::endl;
        return;
    }

    int nx = materialMap.size();
    int ny = materialMap[0].size();

    // VTK 文件头（Structured Grid）
    fout << "# vtk DataFile Version 3.0\n";
    fout << "Material Map from FDTD Simulation\n";
    fout << "ASCII\n";
    fout << "DATASET STRUCTURED_GRID\n";
    fout << "DIMENSIONS " << nx << " " << ny << " 1\n";
    fout << "POINTS " << nx * ny << " float\n";

    // 写入网格点坐标（假设每个格点大小为 dx）
    const double dx = 1.0; // 与你的仿真一致
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
        fout << i * dx << " " << j * dx << " 0.0\n";
        }
    }


    // 写入材料类型数据
    fout << "POINT_DATA " << nx * ny << "\n";
    fout << "SCALARS material_type int 1\n";
    fout << "LOOKUP_TABLE default\n";

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            // 将 MaterialType 转换为整数
            int materialValue = 0;
            if (materialMap[i][j] == MaterialType::VACUUM) {
                materialValue = 0;
            } 
            else if (materialMap[i][j] == MaterialType::DIELECTRIC) {
            materialValue = 1;
            } 
            else if (materialMap[i][j] == MaterialType::METAL) {
            materialValue = 2;
            } 
            else if (materialMap[i][j] == MaterialType::PML) {
            materialValue = 3;
            }
            fout << materialValue << "\n";
        }
    }

    fout.close();
    std::cout << "Saved material map to " << filename << std::endl;
    } 