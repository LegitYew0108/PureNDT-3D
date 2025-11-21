#include "ndt_calculator.hpp"

#include "ndt_core.hpp"
#include "voxel.hpp"
#include <Eigen/Dense>
#include <Eigen/src/Core/util/Constants.h>

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

    // Avoid singular covariance matrix(eq. 6.11)[Magnusson, 2009]
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(
        voxel_pair.second.covariance);
    if (eigensolver.info() != Eigen::Success) {
      // if eigensolver failed
      voxel_pair.second.has_valid_covariance = false;
      continue;
    }

    Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
    Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();

    double lambda_max = eigenvalues(2);

    if (lambda_max < 1e-6) {
      // if max eigenvalue is too small (= points are too dense), skip this
      // voxel.
      voxel_pair.second.has_valid_covariance = false;
      continue;
    }

    double min_lambda = lambda_max / 100.0;
    bool modified = false;

    if (eigenvalues(0) < min_lambda) {
      eigenvalues(0) = min_lambda;
      modified = true;
    }
    if (eigenvalues(1) < min_lambda) {
      eigenvalues(1) = min_lambda;
      modified = true;
    }

    if (modified) {
      Eigen::DiagonalMatrix<double, 3> diagonal_lambdas(
          eigenvalues(0), eigenvalues(1), eigenvalues(2));
      voxel_pair.second.covariance =
          eigenvectors * diagonal_lambdas * eigenvectors.inverse();
    }

    // Perform calculations using constants for numerical stability.(derived
    // from ndt_omp)
    double det = voxel_pair.second.covariance.determinant();
    if (det <= 1e-9) {
      voxel_pair.second.has_valid_covariance = false;
      continue;
    }
    double c_1 = (1.0 - outlier_ratio) / sqrt(pow(2 * M_PI, 3) * det);
    double c_2 = outlier_ratio / pow(resolution_m, 3);
    double d_3 = -std::log(c_2);
    voxel_pair.second.d_1 = -std::log(c_1 + c_2) - d_3;
    voxel_pair.second.d_2 =
        -2 * std::log((-std::log(c_1 * exp(-0.5) + c_2) - d_3) /
                      voxel_pair.second.d_1);

    voxel_pair.second.has_valid_covariance = true;
    log(logger, LogLevel::Debug,
        "[NDTCalculator] Voxel(%d, %d, %d): avg:(%lf, %lf, %lf)",
        voxel_pair.first.x, voxel_pair.first.y, voxel_pair.first.z,
        voxel_pair.second.average.x(), voxel_pair.second.average.y(),
        voxel_pair.second.average.z());
    log(logger, LogLevel::Debug,
        "[NDTCalculator] covariance_mat:\n\t\t(%lf, %lf, %lf)\n\t\t(%lf, %lf, "
        "%lf)\n\t\t(%lf, %lf, %lf)",
        voxel_pair.second.covariance(0), voxel_pair.second.covariance(1),
        voxel_pair.second.covariance(2), voxel_pair.second.covariance(3),
        voxel_pair.second.covariance(4), voxel_pair.second.covariance(5),
        voxel_pair.second.covariance(6), voxel_pair.second.covariance(7),
        voxel_pair.second.covariance(8));
    log(logger, LogLevel::Debug, "[NDTCalculator] d_1: %lf\t\td_2: %lf",
        voxel_pair.second.d_1, voxel_pair.second.d_2);
  }
  log(logger, LogLevel::Info,
      "[NDTCalculator] Completed to calculate NDT. total voxels:%d",
      voxel_grid.get_voxels()->size());
}

} // namespace NDTCalculator
} // namespace PureNDT3D
