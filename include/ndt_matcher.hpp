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
using TransformVec6D = Eigen::Vector<double, 6>;

struct TransformWithVoxelGrid {
  Transform4D transform;
  std::unique_ptr<VoxelGrid> voxel_grid;
};

class NDTMatcher {
public:
  explicit NDTMatcher(const NDTConfig &config);

  TransformType align(const std::vector<Point3D> &source_points,
                      const VoxelGrid &voxel_grid,
                      const TransformType initial_transform);
  const NDTConfig &config_;
  std::unique_ptr<NDTOptimizer> optimizer_;
};
} // namespace PureNDT3D
