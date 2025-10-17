#include "ndt_matcher.hpp"
#include "ndt_core.hpp"
#include "ndt_optimizer.hpp"

namespace PureNDT3D {
Transform4D NDTMatcher::align(const std::vector<Point3D> &source_points,
                              VoxelGrid &voxel_grid,
                              const Transform4D initial_transform) {
  int iteration = 0;
  Transform4D current_transform = initial_transform;
  while (true) {
    TransformWithScore transform_with_score =
        optimizer_.calc_update(source_points, voxel_grid, current_transform);
    current_transform = transform_with_score.transform;
    double score = transform_with_score.score;
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
