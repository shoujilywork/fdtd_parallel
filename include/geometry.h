#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include "utils.h"
// 前向声明 MaterialType 枚举类型

enum MaterialType;

// 创建 V 型天线结构
// 参数说明：
// - centerX, centerY: 天线中心坐标
// - armLength: 每个臂的长度（单位：网格点数）
// - angleDeg: V 型夹角（单位：度）
// - width: 天线宽度（单位：网格点数）
// - leftFacing: 是否开口朝左（true: 左开口，false: 右开口）
void createVShapedAntenna(int centerX, int centerY, int armLength,
                          int angleDeg, int width, bool leftFacing,
                          std::vector<std::vector<MaterialType>>& materialMap);

#endif // GEOMETRY_H