#include "voxel.hpp"
#include "config.hpp"
#include "ndt_calculator.hpp"
#include <unordered_map>

namespace PureNDT3D {
std::unordered_map<VoxelIndex, Voxel, VoxelHash> *VoxelGrid::get_voxels() {
  return &voxels_;
}

Voxel::Voxel()
    : average(Point3D::Zero()), covariance(Covariance3D::Zero()), d_1(0.0),
      d_2(0.0), has_valid_covariance(false) {}

const Voxel *VoxelGrid::get_voxel_const(const VoxelIndex &index) const {
  auto it = voxels_.find(index);
  if (it == voxels_.end()) {
    return nullptr;
  } else {
    return &it->second;
  }
}

const Voxel *VoxelGrid::get_voxel_const(const Point3D &point) const {
  return get_voxel_const(get_point_index(point));
}

VoxelGrid::VoxelGrid(const NDTConfig &config) : config_(config) { ; }

void VoxelGrid::add_points(std::vector<Point3D> points) {
  for (Point3D point : points) {
    VoxelIndex index = get_point_index(point);
    voxels_[index].points.push_back(point);
  }

  // 平均と共分散行列を求める
  NDTCalculator::calc_statistics(config_.logger_, *this, config_.outlier_ratio_,
                                 config_.voxel_resolution_m_);
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

std::vector<Point3D> VoxelGrid::get_average_points() const {
  std::vector<Point3D> average_points;
  for (std::pair<VoxelIndex, Voxel> voxel : voxels_) {
    average_points.push_back(voxel.second.average);
  }
  return average_points;
}
} // namespace PureNDT3D
