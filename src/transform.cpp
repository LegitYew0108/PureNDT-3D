#include "transform.hpp"
namespace {
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
} // namespace

namespace PureNDT3D {
void TransformType::sync_vec_from_mat() {
  vec_.head<3>() = mat_.block<3, 1>(0, 3);
  vec_.tail<3>() = mat_.block<3, 3>(0, 0).eulerAngles(0, 1, 2);
}

void TransformType::sync_mat_from_vec() {
  // vec_ contains [x, y, z, roll, pitch, yaw]
  // We need to create a X-Y-Z rotation matrix
  RotationMatrix3D R_x = get_rotation_mat_x(vec_(3)); // roll
  RotationMatrix3D R_y = get_rotation_mat_y(vec_(4)); // pitch
  RotationMatrix3D R_z = get_rotation_mat_z(vec_(5)); // yaw
  mat_ = Transform4D::Identity();
  mat_.block<3, 3>(0, 0) = R_x * R_y * R_z;
  mat_.block<3, 1>(0, 3) = vec_.head<3>();
}

Transform4D TransformType::se3_exp(const Eigen::Vector<double, 6> &xi) {
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

Eigen::Matrix3d TransformType::skew(const Point3D &p) {
  Eigen::Matrix3d S;
  S << 0, -p.z(), p.y(), p.z(), 0, -p.x(), -p.y(), p.x(), 0;
  return S;
}

TransformType::TransformType(const TransformVec6D &vec) : vec_(vec) {
  sync_mat_from_vec(); // vec から mat を計算
}
TransformType::TransformType(const Transform4D &mat) : mat_(mat) {
  sync_vec_from_mat(); // mat から vec を計算
}

void TransformType::update(TransformVec6D &transform_inc,
                           const double &trans_alpha, const double &rot_alpha) {
  transform_inc.normalize();
  Eigen::Vector<double, 6> alpha_vec;
  alpha_vec << trans_alpha, trans_alpha, trans_alpha, rot_alpha, rot_alpha,
      rot_alpha;
  transform_inc = transform_inc.array() * alpha_vec.array();
  vec_ = vec_ + transform_inc;

  sync_mat_from_vec();
}

// 4x4行列を返す (se3_exp での更新や、最終結果として使う)
const Transform4D &TransformType::get_matrix() const { return mat_; }

// 6Dベクトルを返す (ヤコビアン計算に使う)
const TransformVec6D &TransformType::get_vector() const { return vec_; }

RotationMatrixes TransformType::get_rotation_mat() const {
  RotationMatrixes rot_mats;
  rot_mats.x = get_rotation_mat_x(vec_(3));
  rot_mats.y = get_rotation_mat_y(vec_(4));
  rot_mats.z = get_rotation_mat_z(vec_(5));
  return rot_mats;
}

Point3D TransformType::transform(Eigen::Vector3d point) const {
  Eigen::Vector<double, 4> homogeneous_vec = Eigen::Vector<double, 4>::Ones();
  homogeneous_vec.head<3>() = point;
  return (mat_ * homogeneous_vec).block<3, 1>(0, 0);
}
} // namespace PureNDT3D
