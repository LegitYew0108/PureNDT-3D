#include "ndt_core.hpp"
#include "ndt_matcher.hpp"
#include <cstdarg>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>

namespace PureNDT3D {
NDTCore::NDTCore(const NDTConfig &configs) {
  set_configurations(configs);
  voxel_grid_ = std::make_unique<VoxelGrid>(configs_);
  matcher_ = std::make_unique<NDTMatcher>(configs_);
}
NDTCore::~NDTCore() {}

void log(LoggerCallback logger, LogLevel level, const char *format, ...) {
  if (!logger) {
    return;
  }

  va_list args1;
  va_start(args1, format);

  va_list args2;
  // copy args1 because vsnprintf consume va_list.
  va_copy(args2, args1);
  // get argument size
  int size = vsnprintf(nullptr, 0, format, args2);
  va_end(args2);

  if (size < 0) {
    va_end(args1);
    return;
  }

  // create buffer
  std::vector<char> buffer(static_cast<size_t>(size) + 1);
  vsnprintf(buffer.data(), buffer.size(), format, args1);
  va_end(args1);

  logger(level, std::string(buffer.data(), static_cast<size_t>(size)));
}

void NDTCore::set_configurations(const NDTConfig &configs) {
  configs_ = configs;
}

void NDTCore::check_config(NDTConfig &configs) {
  if (configs.score_threshold_ < 0.0f) {
    throw std::invalid_argument(
        "The score_threshold cannot be set to a negative value.");
  }
  // TODO
}

void NDTCore::add_target_points(const std::vector<Point3D> &points) {
  this->voxel_grid_->add_points(points);
}

void NDTCore::remove_all_target_points() { voxel_grid_->remove_all_points(); }

void NDTCore::replace_target_points(const std::vector<Point3D> &points) {
  remove_all_target_points();
  add_target_points(points);
}

TransformType NDTCore::align(const std::vector<Point3D> &points,
                             const TransformType &initial_transform) {
  return matcher_->align(points, *voxel_grid_, initial_transform);
}
} // namespace PureNDT3D
