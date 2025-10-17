#include "ndt_core.hpp"

namespace PureNDT3D {
NDTCore::NDTCore() {}
NDTCore::NDTCore(const NDTConfig &configs) { set_configurations(configs); }
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

void NDTCore::add_target_points(const std::vector<Point3D> &points) {}

void NDTCore::remove_target_points() {}

void NDTCore::replace_target_points(const std::vector<Point3D> &points) {}

Transform3D NDTCore::align(const std::vector<Point3D> &points,
                           const Transform3D &initial_transform) {}

} // namespace PureNDT3D
