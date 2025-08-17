cmake_minimum_required(VERSION 3.18)
project(fdtd_spp_simulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加可执行文件
add_executable(fdtd_sim
    src/main.cpp
    src/geometry.cpp
    src/materials.cpp
    src/sources.cpp
    src/utils.cpp
    src/pml.cpp
)

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)