#pragma once

#include "config.hpp"
#include "ndt_optimizer.hpp"
#include "voxel.hpp"
#include <eigen3/Eigen/Core>
#include <vector>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform4D = Eigen::Matrix4d;

class NDTMatcher {
public:
  explicit NDTMatcher();

  Transform4D align(const std::vector<Point3D> &source_points,
                    VoxelGrid &voxel_grid, const Transform4D initial_transform);
  NDTConfig config_;
  NDTOptimizer optimizer_;
};
} // namespace PureNDT3D
