#pragma once

#include "config.hpp"
#include "voxel.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform3D = Eigen::Matrix3d;

struct TransformWithScore {
  Transform3D transform;
  double score;
};

class NDTOptimizer {
public:
  explicit NDTOptimizer(const NDTConfig &config);

  Transform3D calc_update(const std::vector<Point3D> &source_points,
                          const VoxelGrid &voxel_grid,
                          const Transform3D &current_transform);

  double get_score(const std::vector<Point3D> &points, VoxelGrid &voxel_grid);

  double ndt_pdf(const Point3D &point, const Voxel &voxel);

  VoxelIndex get_point_index(const Point3D &point);

  NDTConfig config_;
};
} // namespace PureNDT3D
