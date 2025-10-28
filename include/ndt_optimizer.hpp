#pragma once

#include "config.hpp"
#include "voxel.hpp"
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <memory>

namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Transform4D = Eigen::Matrix4d;
using RotationMatrix3D = Eigen::Matrix3d;
using TransformVec6D = Eigen::Vector<double, 6>;

struct RotationMatrixes {
  RotationMatrix3D x;
  RotationMatrix3D y;
  RotationMatrix3D z;
};

Point3D transform(Eigen::Vector3d point, Transform4D transform_mat,
                  RotationMatrixes rotation_mat);

Point3D transform(Eigen::Vector3d point, Transform4D transform_mat);
RotationMatrixes get_rotation_mat(TransformVec6D transform_vec);
RotationMatrix3D get_rotation_mat_x(double theta);
RotationMatrix3D get_rotation_mat_y(double theta);
RotationMatrix3D get_rotation_mat_z(double theta);

struct TransformUpdateType {
  Transform4D transform;
  double score;
};

class NDTOptimizer {
public:
  explicit NDTOptimizer(const NDTConfig &config);

  TransformUpdateType calc_update(const std::vector<Point3D> &source_points,
                                  const VoxelGrid &voxel_grid,
                                  const TransformVec6D &current_transform);

  double get_score(const Point3D &transformed_point, const Voxel &voxel);

protected:
  Eigen::Matrix<double, 3, 6> get_jacobian(const Point3D &point,
                                           const TransformVec6D &rotation_mat);

  Eigen::Matrix<Eigen::Vector3d, 6, 6>
  get_second_order_derivative(const Point3D &point,
                              const TransformVec6D &transform);

  Eigen::Matrix<double, 6, 6> get_hessian(
      const Point3D &points, const double &before_score,
      const Eigen::Matrix<double, 3, 6> &J,
      const Eigen::Matrix<Eigen::Vector3d, 6, 6> &second_order_derivative,
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
