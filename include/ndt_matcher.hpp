#pragma once

#include "config.hpp"
#include "ndt_optimizer.hpp"
#include "voxel.hpp"
#include <eigen3/Eigen/Core>
#include <memory>
#include <vector>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform4D = Eigen::Matrix4d;

class NDTMatcher {
public:
  explicit NDTMatcher(const NDTConfig &config);

  Transform4D align(const std::vector<Point3D> &source_points,
                    std::unique_ptr<VoxelGrid> voxel_grid,
                    const Transform4D initial_transform);
  NDTConfig config_;
  std::unique_ptr<NDTOptimizer> optimizer_;
};
} // namespace PureNDT3D
