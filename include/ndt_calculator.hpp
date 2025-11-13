#pragma once

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "voxel.hpp"
namespace PureNDT3D {
namespace NDTCalculator {
using Point3D = Eigen::Vector3d;
using Covariance3D = Eigen::Matrix3d;
using Transform4D = Eigen::Matrix4d;

void calc_statistics(LoggerCallback logger, VoxelGrid &voxel_grid,
                     double outlier_ratio, double resolution_m);
} // namespace NDTCalculator
} // namespace PureNDT3D
