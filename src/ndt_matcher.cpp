#include "ndt_matcher.hpp"
#include "ndt_core.hpp"
#include "ndt_optimizer.hpp"
#include "voxel.hpp"
#include <memory>

namespace PureNDT3D {
NDTMatcher::NDTMatcher(const NDTConfig &config) : config_(config) {
  optimizer_ = std::make_unique<NDTOptimizer>(config);
}

TransformType NDTMatcher::align(const std::vector<Point3D> &source_points,
                                const VoxelGrid &voxel_grid,
                                const TransformType initial_transform) {
  int iteration = 0;
  TransformType current_transform = initial_transform;

  while (true) {
    double score =
        optimizer_->calc_update(source_points, voxel_grid, current_transform);
    iteration++;
    if (score < config_.score_threshold_) {
      break;
    }
    if (iteration > config_.max_iterations_) {
      break;
    }
  }
  return current_transform;
}
} // namespace PureNDT3D
