#include "ndt_calculator.hpp"
#include "ndt_core.hpp"
#include "voxel.hpp"

namespace PureNDT3D {
namespace NDTCalculator {

void calc_statistics(LoggerCallback logger, VoxelGrid &voxel_grid,
                     double outlier_ratio, double resolution_m) {
  log(logger, LogLevel::Info,
      "[NDTCalculator] Calculate Voxels' Average and Covariance Matrix.");
  // each voxel in voxelgrid
  for (auto &voxel_pair : *voxel_grid.get_voxels()) {
    if (voxel_pair.second.points.empty()) {
      voxel_pair.second.has_valid_covariance = false;
      continue;
    }
    if (voxel_pair.second.points.size() < 5) {
      voxel_pair.second.has_valid_covariance = false;
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
    double d_3 = -std::log(c_2);
    voxel_pair.second.d_1 = -std::log(c_1 + c_2) - d_3;
    voxel_pair.second.d_2 =
        -2 * std::log((-std::log(c_1 * exp(-0.5) + c_2) - d_3) /
                      voxel_pair.second.d_1);

    voxel_pair.second.has_valid_covariance = true;
  }
  log(logger, LogLevel::Info,
      "[NDTCalculator] Completed to calculate NDT. total voxels:%lf",
      voxel_grid.get_voxels()->size());
}

} // namespace NDTCalculator
} // namespace PureNDT3D
