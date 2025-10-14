#include "ndt_core.hpp"

namespace PureNDT3D {
NDTCore::NDTCore() {}
NDTCore::NDTCore(const NDTConfig &configs) {}
NDTCore::~NDTCore() {}

void NDTCore::set_configurations(const NDTConfig &configs) {}

void NDTCore::add_target_points(const std::vector<Point3D> &points) {}

void NDTCore::remove_target_points() {}

void NDTCore::replace_target_points(const std::vector<Point3D> &points) {}

Transform3D NDTCore::align(const std::vector<Point3D> &points,
                           const Transform3D &initial_transform) {}

} // namespace PureNDT3D
