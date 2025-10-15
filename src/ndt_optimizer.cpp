#include "ndt_optimizer.hpp"
#include "voxel.hpp"
#include <cmath>

namespace PureNDT3D {
Transform3D NDTOptimizer::calc_update(const std::vector<Point3D> &source_points,
                                      const VoxelGrid &voxel_grid,
                                      const Transform3D &current_transform) {}

double NDTOptimizer::get_score(const std::vector<Point3D> &points,
                               VoxelGrid &voxel_grid) {
  double score;
  for (Point3D point : points) {
    VoxelIndex index = get_point_index(point);
    // TODO:at関数のtry_catch実装
    Voxel voxel = voxel_grid.get_voxels()->at(index);
    score += ndt_pdf(point, voxel);
  }
  score = -1 * score;
  return score;
}

double NDTOptimizer::ndt_pdf(const Point3D &point, const Voxel &voxel) {
  double prob_density;

  // TODO: inverseの高速化
  prob_density = (point - voxel.average).transpose() *
                 voxel.covariance.inverse() * (point - voxel.average);
  prob_density = std::exp(-1 * prob_density / 2);
  prob_density /= sqrt(pow(M_PI, 3) * voxel.covariance.determinant());
  return prob_density;
}

VoxelIndex NDTOptimizer::get_point_index(const Point3D &point) {
  VoxelIndex index;
  index.x =
      static_cast<int>(std::floor(point(0) / config_.voxel_resolution_m_));
  index.y =
      static_cast<int>(std::floor(point(1) / config_.voxel_resolution_m_));
  index.z =
      static_cast<int>(std::floor(point(2) / config_.voxel_resolution_m_));
  return index;
}

} // namespace PureNDT3D
