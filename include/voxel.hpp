/**
 * @file	voxel.hpp
 * @brief	declalation of Voxel and VoxelGrid.
 * @date	2025/10/10
 * @author	LegitYew0108(Wada Haruto)
 * @note	hash combine funcs are defined here.
 */

#pragma once

#include <eigen3/Eigen/Core>
#include <functional>
#include <unordered_map>

#include <vector>

namespace PureNDT3D {

/**
 * @struct VoxelIndex
 * @brief Index corresponding to each coordinate for using unordered_map in
 * VoxelGrid.
 */
struct VoxelIndex {
  int x, y, z;

  bool operator==(const VoxelIndex &other) const {
    return x == other.x && y == other.y && z == other.z;
  }
};

/**
 * @struct VoxelHash
 * @brief Structure for calculating the hash value of VoxelIndex.
 */
struct VoxelHash {
  // from Boost hash_combine func(64bit ver)
  size_t hash_combine(size_t seed, size_t new_hash) const {
    return new_hash + 0x9e3779b97f4a7c15LLU + (seed << 12) + (seed >> 4);
  }

  size_t operator()(const PureNDT3D::VoxelIndex &i) const {
    size_t seed = std::hash<int>()(i.x);
    seed ^= hash_combine(seed, std::hash<int>()(i.y));
    seed ^= hash_combine(seed, std::hash<int>()(i.z));
    return seed;
  }
};

using Point3D = Eigen::Vector3d;
using Covariance3D = Eigen::Matrix4d;

struct Voxel {
  // for NDT
  Point3D average;
  Covariance3D covariance;

  std::vector<Point3D> points;
};

class VoxelGrid {
public:
  VoxelGrid();

  void add_points(std::vector<Point3D>);

protected:
  std::unordered_map<VoxelIndex, Voxel, VoxelHash> voxels_;
};

} // namespace PureNDT3D
