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

  log(config_.logger_, LogLevel::Info,
      "[NDTMatcher] NDT Matching Loop Start==============================");
  while (true) {
    double score =
        optimizer_->calc_update(source_points, voxel_grid, current_transform);
    log(config_.logger_, LogLevel::Info, "[NDTMatcher] Score: %lf", score);
    iteration++;
    if (score < config_.score_threshold_) {
      log(config_.logger_, LogLevel::Info,
          "[NDTMatcher] Score below the threshold has been detected. finishing "
          "Matching Loop...");
      break;
    }
    if (iteration > config_.max_iterations_) {
      log(config_.logger_, LogLevel::Info,
          "[NDTMatcher] Matching loop iteration was reached max_iterations. "
          "finishing Matching Loop....");
      break;
    }
  }
  return current_transform;
}
} // namespace PureNDT3D
