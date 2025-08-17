// include/types.h
#ifndef TYPES_H
#define TYPES_H

// 材料类型
enum class MaterialType {
    VACUUM,
    DIELECTRIC,
    METAL,
    PML
};

// Drude 模型参数结构体
struct DrudeModel {
    double omega_p;     // 等离子体频率 (rad/s)
    double gamma;       // 碰撞率 (rad/s)
    double eps_inf;     // 高频极限介电常数

    // 默认构造函数
    DrudeModel() : omega_p(0.0), gamma(0.0), eps_inf(1.0) {}

    // 自定义构造函数
    DrudeModel(double wp, double g, double inf)
        : omega_p(wp), gamma(g), eps_inf(inf) {}
};

#endif // TYPES_H