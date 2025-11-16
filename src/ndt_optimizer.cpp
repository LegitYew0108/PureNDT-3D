#include "ndt_optimizer.hpp"

#include <eigen3/Eigen/src/Core/Matrix.h>
#include <omp.h>

#include <cmath>

#include "config.hpp"
#include "ndt_core.hpp"
#include "voxel.hpp"

namespace PureNDT3D {
NDTOptimizer::NDTOptimizer(const NDTConfig &config) : config_(config) { ; }

double NDTOptimizer::calc_update(const std::vector<Point3D> &source_points,
                                 const VoxelGrid &voxel_grid,
                                 TransformType &current_transform) {
  Eigen::Vector<double, 6> gradient = Eigen::Vector<double, 6>::Zero();
  Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
  int valid_point_num = 0;
  double total_score = 0.0;

#pragma omp declare reduction(+ : Eigen::Matrix<double, 6, 6> : omp_out =      \
                                  omp_out + omp_in)                            \
    initializer(omp_priv = omp_orig)
#pragma omp declare reduction(+ : Eigen::Vector<double, 6> : omp_out =         \
                                  omp_out + omp_in)                            \
    initializer(omp_priv = omp_orig)
#pragma omp parallel for reduction(+ : gradient, H, total_score,               \
                                       valid_point_num)
  for (Point3D point : source_points) {
    // 座標変換
    Point3D transformed_point = current_transform.transform(point);

    const Voxel *voxel = voxel_grid.get_voxel_const(transformed_point);
    if (voxel == nullptr) {
      continue;
    }
    if (!voxel->has_valid_covariance) {
      continue;
    }

    double score = get_score(transformed_point, *voxel);

    // その点のヤコビアンを求める
    Eigen::Matrix<double, 3, 6> J = get_jacobian(point, current_transform);

    // 勾配を求めて、総和に加える
    gradient += get_gradient(transformed_point, score, J, *voxel);

    if (config_.use_second_order_derivative_) {
      // 2次の微分項を使う
      Eigen::Matrix<Eigen::Vector3d, 6, 6> second_order_derivative =
          get_second_order_derivative(point, current_transform);
      // ヘッセ行列を求めて、総和に加える
      H += get_hessian(transformed_point, score, J, second_order_derivative,
                       *voxel);
    } else {
      // 2次の微分項を使わない
      // ヘッセ行列を求めて、総和に加える
      H += get_no_second_order_derivative_hessian(transformed_point, score, J,
                                                  *voxel);
    }

    total_score += score;
    valid_point_num++;
  }

  if (valid_point_num == 0) {
    log(config_.logger_, LogLevel::Error,
        "No valid point detected. ignoring result...");
    return 0.0;
  }

  TransformVec6D transform_inc =
      -(H + config_.levenberg_marquardt_lambda_ *
                Eigen::Matrix<double, 6, 6>::Identity())
           .ldlt()
           .solve(gradient);

  current_transform.update(transform_inc);

  return total_score / valid_point_num;
}

double NDTOptimizer::get_score(const Point3D &point, const Voxel &voxel) {
  Eigen::Vector3d diff = point - voxel.average;
  double exponent =
      -0.5 * voxel.d_2 *
      (diff.transpose() * voxel.covariance.inverse() * diff)(0, 0);
  double score = -voxel.d_1 * exp(exponent);
  return score;
}

Eigen::Matrix<double, 3, 6>
NDTOptimizer::get_jacobian(const Point3D &point,
                           const TransformType &transform) {
  // Get rotation angles from the transform vector, which is [x, y, z, yaw, pitch, roll]
  // vec[3] -> yaw (rot_z)
  // vec[4] -> pitch (rot_y)
  // vec[5] -> roll (rot_x)
  double rot_z = transform.get_vector()[3]; // yaw
  double rot_y = transform.get_vector()[4]; // pitch
  double rot_x = transform.get_vector()[5]; // roll

  // Pre-calculate sin and cos
  double s_x = sin(rot_x);
  double c_x = cos(rot_x);
  double s_y = sin(rot_y);
  double c_y = cos(rot_y);
  double s_z = sin(rot_z);
  double c_z = cos(rot_z);

  // Rotation matrices (for Z-Y-X Euler angles)
  Eigen::Matrix3d R_x, R_y, R_z;
  R_x << 1, 0, 0, 0, c_x, -s_x, 0, s_x, c_x;
  R_y << c_y, 0, s_y, 0, 1, 0, -s_y, 0, c_y;
  R_z << c_z, -s_z, 0, s_z, c_z, 0, 0, 0, 1;

  // Derivative of rotation matrices
  Eigen::Matrix3d dR_dx, dR_dy, dR_dz;
  dR_dx << 0, 0, 0, 0, -s_x, -c_x, 0, c_x, -s_x; // Derivative w.r.t. roll
  dR_dy << -s_y, 0, c_y, 0, 0, 0, -c_y, 0, -s_y; // Derivative w.r.t. pitch
  dR_dz << -s_z, -c_z, 0, c_z, -s_z, 0, 0, 0, 0; // Derivative w.r.t. yaw

  // Jacobian calculation using matrix chain rule for a Z-Y-X rotation
  // d(R*p)/d(roll)
  Eigen::Vector3d J_rot_x = R_z * R_y * dR_dx * point;
  // d(R*p)/d(pitch)
  Eigen::Vector3d J_rot_y = R_z * dR_dy * R_x * point;
  // d(R*p)/d(yaw)
  Eigen::Vector3d J_rot_z = dR_dz * R_y * R_x * point;

  Eigen::Matrix<double, 3, 6> Jacobian;
  Jacobian.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  // The columns of the Jacobian must match the order in the transform vector: [x, y, z, yaw, pitch, roll]
  Jacobian.col(3) = J_rot_z; // Partial derivative with respect to yaw
  Jacobian.col(4) = J_rot_y; // Partial derivative with respect to pitch
  Jacobian.col(5) = J_rot_x; // Partial derivative with respect to roll

  return Jacobian;
}

Eigen::Matrix<Eigen::Vector3d, 6, 6>
NDTOptimizer::get_second_order_derivative(const Point3D &point,
                                          const TransformType &transform) {
  // Correctly assign Euler angles: [x, y, z, yaw, pitch, roll]
  double rot_z = transform.get_vector()[3]; // yaw
  double rot_y = transform.get_vector()[4]; // pitch
  double rot_x = transform.get_vector()[5]; // roll
  double s_x = sin(rot_x);
  double c_x = cos(rot_x);
  double s_y = sin(rot_y);
  double c_y = cos(rot_y);
  double s_z = sin(rot_z);
  double c_z = cos(rot_z);
  Eigen::Vector3d a, b, c, d, e, f;
  a << 0,
      point[0] * (-c_x * s_z - s_x * s_y * c_z) +
          point[1] * (-c_x * c_z + s_x * s_y * s_z) + point[2] * (s_x * c_y),
      point[0] * (-s_x * s_z + c_x * s_y * c_z) +
          point[1] * (-c_x * s_y * s_z - s_x * c_z) + point[2] * (-c_x * c_y);

  b << 0,
      point[0] * (c_x * c_y * c_z) + point[1] * (-c_x * c_y * s_z) +
          point[2] * (c_x * s_y),
      point[0] * (s_x * c_y * c_z) + point[1] * (-s_x * c_y * s_z) +
          point[2] * (s_x * s_y);

  c << 0,
      point[0] * (-s_x * c_z - c_x * s_y * s_z) +
          point[1] * (-s_x * s_z - c_x * s_y * c_z),
      point[0] * (c_x * c_z - s_x * s_y * s_z) +
          point[1] * (-s_x * s_y * c_z - c_x * s_z);

  d << point[0] * (-c_y * c_z) + point[1] * (c_y * s_z) + point[2] * (-s_y),
      point[0] * (-s_x * s_y * c_z) + point[1] * (s_x * s_y * s_z) +
          point[2] * (s_x * c_y),
      point[0] * (c_x * s_y * c_z) + point[1] * (-c_x * s_y * s_z) +
          point[2] * (-c_x * c_y);

  e << point[0] * (s_y * s_z) + point[1] * (s_y * c_z),
      point[0] * (-s_x * c_y * s_z) + point[1] * (-s_x * c_y * c_z),
      point[0] * (c_x * c_y * s_z) + point[1] * (c_x * c_y * c_z);

  f << point[0] * (-c_y * c_z) + point[1] * (c_y * s_z),
      point[0] * (-c_x * s_z - s_x * s_y * c_z) +
          point[1] * (-c_x * c_z + s_x * s_y * s_z),
      point[0] * (-s_x * s_z + c_x * s_y * c_z) +
          point[1] * (-c_x * s_y * s_z - s_x * c_z);

  Eigen::Matrix<Eigen::Vector3d, 6, 6> H_mat;
  H_mat.setConstant(Eigen::Vector3d::Zero());
  H_mat.block<3, 3>(3, 3) << a, b, c, b, d, e, c, e, f;
  return H_mat;
}

Eigen::Vector<double, 6> NDTOptimizer::get_gradient(
    const Point3D &transformed_point, const double &before_score,
    const Eigen::Matrix<double, 3, 6> &J, const Voxel &voxel) {
  Point3D diff = transformed_point - voxel.average;
  double exponent =
      -0.5 * voxel.d_2 * diff.transpose() * voxel.covariance.inverse() * diff;
  Eigen::Vector<double, 6> gradient = voxel.d_1 * voxel.d_2 * diff.transpose() *
                                      voxel.covariance.inverse() * J *
                                      exp(exponent);
  return gradient;
}

Eigen::Matrix<double, 6, 6>
NDTOptimizer::get_no_second_order_derivative_hessian(
    const Point3D &transformed_point, const double &before_score,
    const Eigen::Matrix<double, 3, 6> &J, const Voxel &voxel) {
  Point3D diff = transformed_point - voxel.average;
  double exponent =
      -0.5 * voxel.d_2 * diff.transpose() * voxel.covariance.inverse() * diff;
  double coeffs = voxel.d_1 * voxel.d_2 * exp(exponent);

  Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();

  for (int i = 0; i < 6; i++) {
    double a =
        diff.transpose() * voxel.covariance.inverse() * J.block<3, 1>(0, i);
    for (int j = 0; j < 6; j++) {
      double b =
          diff.transpose() * voxel.covariance.inverse() * J.block<3, 1>(0, j);
      H(i, j) = coeffs * (-voxel.d_2 * a * b + J.block<3, 1>(0, j).transpose() *
                                                   voxel.covariance.inverse() *
                                                   J.block<3, 1>(0, i));
    }
  }

  return H;
}

Eigen::Matrix<double, 6, 6> NDTOptimizer::get_hessian(
    const Point3D &transformed_point, const double &before_score,
    const Eigen::Matrix<double, 3, 6> &J,
    const Eigen::Matrix<Eigen::Vector3d, 6, 6> &second_order_derivative,
    const Voxel &voxel) {
  Point3D diff = transformed_point - voxel.average;
  double exponent =
      -0.5 * voxel.d_2 * diff.transpose() * voxel.covariance.inverse() * diff;
  double coeffs = voxel.d_1 * voxel.d_2 * exp(exponent);

  Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();

  for (int i = 0; i < 6; i++) {
    double a =
        diff.transpose() * voxel.covariance.inverse() * J.block<3, 1>(0, i);
    for (int j = 0; j < 6; j++) {
      double b =
          diff.transpose() * voxel.covariance.inverse() * J.block<3, 1>(0, j);
      H(i, j) = coeffs * (-voxel.d_2 * a * b +
                          diff.transpose() * voxel.covariance.inverse() *
                              second_order_derivative(i, j) +
                          J.block<3, 1>(0, j).transpose() *
                              voxel.covariance.inverse() * J.block<3, 1>(0, i));
    }
  }

  return H;
}

} // namespace PureNDT3D
