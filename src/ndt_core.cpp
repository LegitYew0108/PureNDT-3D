#include "ndt_core.hpp"
#include <filesystem>
#include <memory>

namespace PureNDT3D {
NDTCore::NDTCore(const NDTConfig &configs) {
  set_configurations(configs);
  voxel_grid_ = std::make_unique<VoxelGrid>(configs);
  matcher_ = std::make_unique<NDTMatcher>(configs);
}
NDTCore::~NDTCore() {}

void log(LoggerCallback logger, LogLevel level, const char *message, ...) {
  if (logger) {
    logger(level, message);
  }
}

void NDTCore::set_configurations(const NDTConfig &configs) {
  configs_ = configs;
}

void check_config(NDTConfig &configs) {
  if (configs.score_threshold_ < 0.0f) {
    throw;
  }
}

void NDTCore::add_target_points(const std::vector<Point3D> &points) {
  this->voxel_grid_->add_points(points);
}

void NDTCore::remove_all_target_points() {}

void NDTCore::replace_target_points(const std::vector<Point3D> &points) {
  remove_all_target_points();
  add_target_points(points);
}

Transform4D NDTCore::align(const std::vector<Point3D> &points,
                           const Transform4D &initial_transform) {
  return matcher_->align(points, std::move(voxel_grid_), initial_transform);
}
} // namespace PureNDT3D
