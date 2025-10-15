#pragma once

#include "voxel.hpp"
namespace PureNDT3D {
namespace NDTCalculator {
using Point3D = Eigen::Vector3d;
using Covariance3D = Eigen::Matrix3d;
using Transform3D = Eigen::Matrix3d;

void calc_statistics(VoxelGrid &voxel_grid);
} // namespace NDTCalculator
} // namespace PureNDT3D
