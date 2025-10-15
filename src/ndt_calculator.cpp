#include "ndt_calculator.hpp"
#include "voxel.hpp"

namespace PureNDT3D {
namespace NDTCalculator {

void calc_statistics(VoxelGrid &voxel_grid) {
  // each voxel in voxelgrid
  for (auto &voxel_pair : *voxel_grid.get_voxels()) {
    if (voxel_pair.second.points.empty()) {
      continue;
    }

    // get points' average vector
    Point3D sum_point = Point3D::Zero();
    // each points in one voxel.
    for (Point3D point : voxel_pair.second.points) {
      sum_point += point;
    }
    voxel_pair.second.average = sum_point / voxel_pair.second.points.size();

    // get points' covariance matrix
    Covariance3D covariance_sum = Covariance3D::Zero();
    // each points in one voxel.
    for (Point3D point : voxel_pair.second.points) {
      Point3D diff_avg = point - voxel_pair.second.average;
      covariance_sum += diff_avg * diff_avg.transpose();
    }
    voxel_pair.second.covariance =
        covariance_sum / voxel_pair.second.points.size();
  }
}
} // namespace NDTCalculator
} // namespace PureNDT3D
