#include "ndt_optimizer.hpp"
#include "voxel.hpp"
#include <cmath>
#include <eigen3/Eigen/src/Core/Matrix.h>

namespace PureNDT3D {

Point3D transform(Eigen::Vector3d point, TransformVec6D transform_vec,
                  RotationMatrixes rotation_mat) {
  return rotation_mat.z * rotation_mat.y * rotation_mat.x * point +
         transform_vec.block<3, 1>(0, 0);
}
Point3D transform(Eigen::Vector3d point, Transform4D transform_mat) {
  Eigen::Vector4d homogeneous_point = Eigen::Vector4d::Ones();
  homogeneous_point.block<3, 1>(0, 0) = point;
  return (transform_mat * homogeneous_point).block<3, 1>(0, 0);
}

RotationMatrixes get_rotation_mat(TransformVec6D transform_vec) {
  RotationMatrixes rot_matrix = RotationMatrixes();
  rot_matrix.x = get_rotation_mat_x(transform_vec[3]);
  rot_matrix.y = get_rotation_mat_y(transform_vec[4]);
  rot_matrix.z = get_rotation_mat_z(transform_vec[5]);
  return rot_matrix;
}
RotationMatrix3D get_rotation_mat_x(double theta) {
  RotationMatrix3D x_mat = RotationMatrix3D::Zero();
  x_mat << 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta);
  return x_mat;
}
RotationMatrix3D get_rotation_mat_y(double theta) {
  RotationMatrix3D y_mat = RotationMatrix3D::Zero();
  y_mat << cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta);
  return y_mat;
}
RotationMatrix3D get_rotation_mat_z(double theta) {
  RotationMatrix3D z_mat = RotationMatrix3D::Zero();
  z_mat << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
  return z_mat;
}

NDTOptimizer::NDTOptimizer(const NDTConfig &config) : config_(config) { ; }

TransformUpdateType
NDTOptimizer::calc_update(const std::vector<Point3D> &source_points,
                          const VoxelGrid &voxel_grid,
                          const TransformVec6D &current_transform_vec) {

  // TODO
  RotationMatrixes R = get_rotation_mat(current_transform_vec);
  Eigen::Vector<double, 6> gradient = Eigen::Vector<double, 6>::Zero();
  Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();

  for (Point3D point : source_points) {
    // 座標変換
    Point3D transformed_point = transform(point, current_transform_vec, R);

    const Voxel *voxel = voxel_grid.get_voxel_const(transformed_point);
    if (voxel == nullptr) {
      continue;
    }

    double score = get_score(transformed_point, *voxel);

    // その点のヤコビアンを求める
    Eigen::Matrix<double, 3, 6> J = get_jacobian(point, current_transform_vec);

    // 勾配を求めて、総和に加える
    gradient += get_gradient(transformed_point, score, J, *voxel);

    if (config_.use_second_order_derivative_) {
      // 2次の微分項を使う
      Eigen::Matrix<Eigen::Vector3d, 6, 6> second_order_derivative =
          get_second_order_derivative(point, current_transform_vec);
      // ヘッセ行列を求めて、総和に加える
      H += get_hessian(transformed_point, score, J, second_order_derivative,
                       *voxel);
    } else {
      // 2次の微分項を使わない
      // ヘッセ行列を求めて、総和に加える
      H += get_no_second_order_derivative_hessian(transformed_point, score, J,
                                                  *voxel);
    }
  }

  TransformVec6D transform_inc =
      -(H + config_.levenberg_marquardt_lambda_ *
                Eigen::Matrix<double, 6, 6>::Identity())
           .ldlt()
           .solve(gradient);

  Transform4D transform_inc_mat = se3_exp(transform_inc);
  Transform4D current_transform_mat = Transform4D::Zero();
  current_transform_mat.block<3, 3>(0, 0) = R.z * R.y * R.x;
  current_transform_mat.block<3, 1>(0, 3) =
      current_transform_vec.block<3, 1>(0, 0);
  current_transform_mat(3, 3) = 1;
  TransformDatas new_transform;
  new_transform.mat = transform_inc_mat * current_transform_mat;
  // TODO
  new_transform.vec.head<3>() = new_transform.mat.block<3, 1>(0, 3);
  Eigen::Vector3d euler =
      new_transform.mat.block<3, 3>(0, 0).eulerAngles(2, 1, 0);
  new_transform.vec.tail<3>() = euler;
  double total_score = 0.0;

  for (Point3D point : source_points) {
    Point3D transformed_point = transform(point, new_transform.mat);

    const Voxel *voxel = voxel_grid.get_voxel_const(transformed_point);
    if (voxel == nullptr) {
      continue;
    }

    total_score += get_score(transformed_point, *voxel);
  }

  TransformUpdateType update = TransformUpdateType{new_transform, total_score};
  return update;
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
                           const TransformVec6D &transform) {
  double rot_x = transform[3];
  double rot_y = transform[4];
  double rot_z = transform[5];
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
                                          const TransformVec6D &transform) {
  double rot_x = transform[3];
  double rot_y = transform[4];
  double rot_z = transform[5];
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

Transform4D NDTOptimizer::se3_exp(const Eigen::Vector<double, 6> &xi) {
  Eigen::Vector3d rho = xi.head<3>();
  Eigen::Vector3d phi = xi.tail<3>();

  double theta = phi.norm();
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d V = Eigen::Matrix3d::Identity();

  if (theta < 1e-5) {
    // 小さい角度のときはテイラー展開による近似
    R = Eigen::Matrix3d::Identity() + skew(phi);
    V = Eigen::Matrix3d::Identity() + 0.5 * skew(phi);
  } else {
    Eigen::Matrix3d phi_hat = skew(phi);
    R = Eigen::Matrix3d::Identity() + (sin(theta) / theta) * phi_hat +
        ((1 - cos(theta)) / (theta * theta)) * (phi_hat * phi_hat);

    V = Eigen::Matrix3d::Identity() +
        ((1 - cos(theta)) / (theta * theta)) * phi_hat +
        ((theta - sin(theta)) / (theta * theta * theta)) * (phi_hat * phi_hat);
  }

  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3, 3>(0, 0) = R;
  T.block<3, 1>(0, 3) = V * rho;

  return T;
}

Eigen::Matrix3d NDTOptimizer::skew(const Point3D &p) {
  Eigen::Matrix3d S;
  S << 0, -p.z(), p.y(), p.z(), 0, -p.x(), -p.y(), p.x(), 0;
  return S;
}

} // namespace PureNDT3D
