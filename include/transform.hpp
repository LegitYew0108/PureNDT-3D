#pragma once
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"

using Point3D = Eigen::Vector3d;
using Transform4D = Eigen::Matrix4d;
using RotationMatrix3D = Eigen::Matrix3d;
using TransformVec6D = Eigen::Vector<double, 6>;

namespace PureNDT3D {
struct RotationMatrixes {
  RotationMatrix3D x;
  RotationMatrix3D y;
  RotationMatrix3D z;
};

// このような Transform クラスを新しく定義する
class TransformType {
private:
  Transform4D mat_;    // 4x4 行列
  TransformVec6D vec_; // 6D ベクトル (tx, ty, tz, rx, ry, rz)

  void sync_vec_from_mat();

  void sync_mat_from_vec();

  static Transform4D se3_exp(const TransformVec6D &xi);

  static Eigen::Matrix3d skew(const Point3D &p);

public:
  explicit TransformType(const TransformVec6D &vec);
  explicit TransformType(const Transform4D &mat);

  void update(TransformVec6D &transform_inc);

  // 4x4行列を返す (se3_exp での更新や、最終結果として使う)
  const Transform4D &get_matrix() const;

  // 6Dベクトルを返す (ヤコビアン計算に使う)
  const TransformVec6D &get_vector() const;

  RotationMatrixes get_rotation_mat() const;

  Point3D transform(Eigen::Vector3d point) const;
};
} // namespace PureNDT3D
