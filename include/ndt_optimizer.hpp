#pragma once

#include "voxel.hpp"
#include <eigen3/Eigen/Core>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform3D = Eigen::Matrix4d;

class NDTOptimizer {
public:
  explicit NDTOptimizer();

  Transform3D calc_update(const std::vector<Point3D> &source_points,
                          const VoxelGrid &voxel_grid,
                          const Transform3D &current_transform);
};
} // namespace PureNDT3D
