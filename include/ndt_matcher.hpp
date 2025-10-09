#pragma once

#include "voxel.hpp"
#include <eigen3/Eigen/Core>
#include <vector>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform3D = Eigen::Matrix4d;

class NDTMatcher {
public:
  Transform3D align(const std::vector<Point3D> &source_points,
                    const VoxelGrid &voxel_grids,
                    const Transform3D initial_transform);
};
} // namespace PureNDT3D
