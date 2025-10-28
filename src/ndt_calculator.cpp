#include "ndt_calculator.hpp"
#include "voxel.hpp"

namespace PureNDT3D {
namespace NDTCalculator {

void calc_statistics(VoxelGrid &voxel_grid, double outlier_ratio,
                     double resolution_m) {
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

    double det = voxel_pair.second.covariance.determinant();
    double c_1 = (1.0 - outlier_ratio) / sqrt(pow(2 * M_PI, 3) * det);
    double c_2 = 1.0 / pow(resolution_m, 3);
    double d_3 = -log(c_2);
    voxel_pair.second.d_1 = -log(c_1 + c_2) - d_3;
    voxel_pair.second.d_2 =
        -2 * log((-log(c_1 * exp(-0.5) + c_2) - d_3) / voxel_pair.second.d_1);
  }
}

} // namespace NDTCalculator
} // namespace PureNDT3D
