#include "ndt_optimizer.hpp"
#include "voxel.hpp"
#include <cmath>

namespace PureNDT3D {

Point3D transform(Eigen::Vector3d point, Transform4D transform_mat) {
  Eigen::Vector4d homogeneous_point = Eigen::Vector4d::Ones();
  homogeneous_point.block<3, 1>(0, 0) = point;
  return (transform_mat * homogeneous_point).block<3, 1>(0, 0);
}

TransformUpdateType
NDTOptimizer::calc_update(const std::vector<Point3D> &source_points,
                          const VoxelGrid &voxel_grid,
                          const Transform4D &current_transform) {
  // 同次座標変換行列から回転行列部を抜き出す
  RotationMatrix3D R = current_transform.block<3, 3>(0, 0);
  Eigen::Vector<double, 6> gradient = Eigen::Vector<double, 6>::Zero();
  Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();

  for (Point3D point : source_points) {
    // 座標変換
    Point3D transformed_point = transform(point, current_transform);

    // TODO:nullptrのときの処理
    const Voxel *voxel = voxel_grid.get_voxel_const(transformed_point);

    double score = get_score(transformed_point, *voxel);

    // その点のヤコビアンを求める
    Eigen::Matrix<double, 3, 6> J = get_jacobian(point, R);

    // 勾配を求めて、総和に加える
    gradient += get_gradient(transformed_point, score, J, *voxel);

    // ヘッセ行列を求めて、総和に加える
    H += get_hessian(transformed_point, score, J, *voxel);
  }

  Eigen::Vector<double, 6> transform_inc = -H.ldlt().solve(gradient);

  Transform4D transform_inc_mat = se3_exp(transform_inc);
  Transform4D new_transform = transform_inc_mat * current_transform;
  double total_score = 0.0;

  for (Point3D point : source_points) {
    Point3D transformed_point = transform(point, new_transform);

    // TODO:nullptrのときの処理
    const Voxel *voxel = voxel_grid.get_voxel_const(transformed_point);

    total_score += get_score(transformed_point, *voxel);
  }

  TransformUpdateType update = TransformUpdateType{new_transform, total_score};
  return update;
}

double NDTOptimizer::get_score(const Point3D &transformed_point,
                               const Voxel &voxel) {
  double score = -ndt_pdf(transformed_point, voxel);
  return score;
}

double NDTOptimizer::ndt_pdf(const Point3D &point, const Voxel &voxel) {
  Eigen::Vector3d diff = point - voxel.average;
  double exponent = -0.5 * diff.transpose() * voxel.covariance.inverse() * diff;
  double det = voxel.covariance.determinant();
  double coeff = 1.0 / sqrt(pow(2 * M_PI, 3) * det);
  return coeff * std::exp(exponent);
}

Eigen::Matrix<double, 3, 6>
NDTOptimizer::get_jacobian(const Point3D &point, const RotationMatrix3D &R) {
  Eigen::Matrix<double, 3, 6> Jacobian;
  Jacobian.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  Jacobian.block<3, 3>(0, 3) = -R * skew(point);
  return Jacobian;
}

Eigen::Vector<double, 6> NDTOptimizer::get_gradient(
    const Point3D &transformed_point, const double &before_score,
    const Eigen::Matrix<double, 3, 6> &J, const Voxel &voxel) {

  return -before_score * J.transpose() * voxel.covariance.inverse() *
         (transformed_point - voxel.average);
}

Eigen::Matrix<double, 6, 6> NDTOptimizer::get_hessian(
    const Point3D &transformed_point, const double &before_score,
    const Eigen::Matrix<double, 3, 6> &J, const Voxel &voxel) {
  Eigen::Vector<double, 6> j_sigma_d_vec = J.transpose() *
                                           voxel.covariance.inverse() *
                                           (transformed_point - voxel.average);
  Eigen::Matrix<double, 6, 6> j_sigma_d_mat =
      j_sigma_d_vec * j_sigma_d_vec.transpose();

  Eigen::Matrix<double, 6, 6> H =
      j_sigma_d_mat - J.transpose() * voxel.covariance.inverse() * J;
  H *= before_score;

  return H;
}

Transform4D NDTOptimizer::se3_exp(const Eigen::Vector<double, 6> &xi) {
  Eigen::Vector3d rho = xi.head<3>();
  Eigen::Vector3d phi = xi.tail<3>();

  double theta = phi.norm();
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d V = Eigen::Matrix3d::Identity();

  if (theta < 1e-8) {
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
