#pragma once

#include "config.hpp"
#include "voxel.hpp"
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform4D = Eigen::Matrix4d;
using RotationMatrix3D = Eigen::Matrix3d;

Point3D transform(Eigen::Vector3d point, Transform4D transform_mat);

class NDTOptimizer {
public:
  explicit NDTOptimizer(const NDTConfig &config);

  Transform4D
  calc_update(const std::vector<std::pair<Point3D, double>> &source_points,
              VoxelGrid &voxel_grid, const Transform4D &current_transform);

  double get_score(const Point3D &transformed_point, Voxel &voxel);

protected:
  double ndt_pdf(const Point3D &point, const Voxel &voxel);

  VoxelIndex get_point_index(const Point3D &point);

  Eigen::Matrix<double, 3, 6>
  get_jacobian(const Point3D &point, const RotationMatrix3D &rotation_mat);

  Eigen::Matrix<double, 6, 6> get_hessian(const Point3D &points,
                                          const double &before_score,
                                          const Eigen::Matrix<double, 3, 6> &J,
                                          const Voxel &voxel);

  Eigen::Vector<double, 6> get_gradient(const Point3D &points,
                                        const double &before_score,
                                        const Eigen::Matrix<double, 3, 6> &J,
                                        const Voxel &voxel);

  Transform4D se3_exp(const Eigen::Vector<double, 6> &xi);

  Eigen::Matrix3d skew(const Point3D &p);

  NDTConfig config_;
};
} // namespace PureNDT3D
