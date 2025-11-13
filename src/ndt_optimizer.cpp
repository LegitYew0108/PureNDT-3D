#include "ndt_optimizer.hpp"
#include "voxel.hpp"
#include <cmath>
#include <eigen3/Eigen/src/Core/Matrix.h>

namespace PureNDT3D {
NDTOptimizer::NDTOptimizer(const NDTConfig &config) : config_(config) { ; }

double NDTOptimizer::calc_update(const std::vector<Point3D> &source_points,
                                 const VoxelGrid &voxel_grid,
                                 TransformType &current_transform) {
  Eigen::Vector<double, 6> gradient = Eigen::Vector<double, 6>::Zero();
  Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
  int valid_point_num = 0;
  double total_score = 0.0;

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
  double rot_x = transform.get_vector()[3];
  double rot_y = transform.get_vector()[4];
  double rot_z = transform.get_vector()[5];
  double s_x = sin(rot_x);
  double c_x = cos(rot_x);
  double s_y = sin(rot_y);
  double c_y = cos(rot_y);
  double s_z = sin(rot_z);
  double c_z = cos(rot_z);
  double a = point[0] * (-s_x * s_z + c_x * s_y * c_z) +
             point[1] * (-s_x * c_z - c_x * s_y * s_z) +
             point[2] * (-c_x * c_y);
  double b = point[0] * (c_x * s_z + s_x * s_y * c_z) +
             point[1] * (-s_x * s_y * s_z + c_x * c_z) +
             point[2] * (-s_x * c_y);
  double c = point[0] * (-s_y * c_z) + point[1] * (s_y * s_z) + point[2] * c_y;
  double d = point[0] * (s_x * c_y * c_z) + point[1] * (-s_x * c_y * s_z) +
             point[2] * (s_x * s_y);
  double e = point[0] * (-c_x * c_y * c_z) + point[1] * (c_x * c_y * s_z) +
             point[2] * -c_x * s_y;
  double f = point[0] * -c_y * s_z + point[1] * -c_y * c_z;
  double g = point[0] * (c_x * c_z - s_x * s_y * s_z) +
             point[1] * (-c_x * s_z - s_x * s_y * c_z);
  double h = point[0] * (s_x * c_z + c_x * s_y * s_z) +
             point[1] * (c_x * s_y * c_z - s_x * s_z);

  Eigen::Matrix<double, 3, 6> Jacobian;
  Jacobian.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  Jacobian.block<3, 3>(0, 3) << 0, c, f, a, d, g, b, e, h;
  return Jacobian;
}

Eigen::Matrix<Eigen::Vector3d, 6, 6>
NDTOptimizer::get_second_order_derivative(const Point3D &point,
                                          const TransformType &transform) {
  double rot_x = transform.get_vector()[3];
  double rot_y = transform.get_vector()[4];
  double rot_z = transform.get_vector()[5];
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
