#pragma once

#include <eigen3/Eigen/Core>
#include <unordered_map>
#include <vector>
namespace PureNDT3D {
using Point3D = Eigen::Vector3d;
using Covariance3D = Eigen::Matrix4d;

struct Voxel {
  Point3D average;
  Covariance3D covariance;
  std::vector<Point3D> points;
};

class VoxelGrid {
public:
  void add_points(std::vector<Point3D>);

protected:
  std::unordered_map<Point3D, Voxel> voxels;
};

} // namespace PureNDT3D
