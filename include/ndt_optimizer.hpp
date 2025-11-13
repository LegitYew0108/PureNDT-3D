#pragma once

#include "config.hpp"
#include "eigen3/Eigen/Cholesky"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"
#include "eigen3/Eigen/LU"
#include "transform.hpp"
#include "voxel.hpp"
#include <memory>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform4D = Eigen::Matrix4d;
using RotationMatrix3D = Eigen::Matrix3d;
using TransformVec6D = Eigen::Vector<double, 6>;

class NDTOptimizer {
public:
  explicit NDTOptimizer(const NDTConfig &config);

  double calc_update(const std::vector<Point3D> &source_points,
                     const VoxelGrid &voxel_grid,
                     TransformType &current_transform);

  double get_score(const Point3D &transformed_point, const Voxel &voxel);

protected:
  Eigen::Matrix<double, 3, 6> get_jacobian(const Point3D &point,
                                           const TransformType &rotation_mat);

  Eigen::Matrix<Eigen::Vector3d, 6, 6>
  get_second_order_derivative(const Point3D &point,
                              const TransformType &transform);

  Eigen::Matrix<double, 6, 6> get_no_second_order_derivative_hessian(
      const Point3D &transformed_point, const double &before_score,
      const Eigen::Matrix<double, 3, 6> &J, const Voxel &voxel);

  Eigen::Matrix<double, 6, 6> get_hessian(
      const Point3D &points, const double &before_score,
      const Eigen::Matrix<double, 3, 6> &J,
      const Eigen::Matrix<Eigen::Vector3d, 6, 6> &second_order_derivative,
      const Voxel &voxel);

  Eigen::Vector<double, 6> get_gradient(const Point3D &points,
                                        const double &before_score,
                                        const Eigen::Matrix<double, 3, 6> &J,
                                        const Voxel &voxel);

  const NDTConfig &config_;
};
} // namespace PureNDT3D
