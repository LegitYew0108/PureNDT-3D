#pragma once

#include <eigen3/Eigen/Core>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform3D = Eigen::Matrix4d;

class NDTOptimizer {
public:
  Transform3D calc_update();
};
} // namespace PureNDT3D
