#include "voxel.hpp"
#include "config.hpp"
#include "ndt_calculator.hpp"
#include <unordered_map>

namespace PureNDT3D {
std::unordered_map<VoxelIndex, Voxel, VoxelHash> *VoxelGrid::get_voxels() {
  return &voxels_;
}

const Voxel VoxelGrid::get_voxel_const(const VoxelIndex &index) const {
  // TODO: try_catch
  const Voxel voxel = voxels_.at(index);
  return voxel;
}

const Voxel VoxelGrid::get_voxel_const(const Point3D &point) const {
  return get_voxel_const(get_point_index(point));
}

VoxelGrid::VoxelGrid(const NDTConfig &config) { config_ = config; }

void VoxelGrid::add_points(std::vector<Point3D> points) {
  for (Point3D point : points) {
    VoxelIndex index = get_point_index(point);
    if (this->voxels_.count(index) == 1) {
      Voxel voxel = voxels_.at(index);
      voxel.points.push_back(point);
    } else {
      // 同じindexのvoxelが存在しなかった時
      // 新しいvoxelを作り、挿入する
      Voxel new_voxel;
      new_voxel.average = Point3D::Zero();
      new_voxel.covariance = Covariance3D::Zero();
      new_voxel.points.push_back(point);
      this->voxels_.emplace(index, new_voxel);
    }
  }

  // 平均と共分散行列を求める
  NDTCalculator::calc_statistics(*this);
}

void VoxelGrid::remove_all_points() { voxels_.clear(); }

VoxelIndex VoxelGrid::get_point_index(const Point3D &point) const {
  VoxelIndex index;
  index.x =
      static_cast<int>(std::floor(point(0) / config_.voxel_resolution_m_));
  index.y =
      static_cast<int>(std::floor(point(1) / config_.voxel_resolution_m_));
  index.z =
      static_cast<int>(std::floor(point(2) / config_.voxel_resolution_m_));
  return index;
}
} // namespace PureNDT3D
