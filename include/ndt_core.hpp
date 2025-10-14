#pragma once

#include "config.hpp"
#include "eigen3/Eigen/Core"
#include "ndt_matcher.hpp"
#include "voxel.hpp"
#include <memory>
#include <vector>

namespace PureNDT3D {

using Point3D = Eigen::Vector3d;
using Covariance3D = Eigen::Matrix4d;
using Transform3D = Eigen::Matrix4d;

class NDTCore {
public:
  /**
   * @brief default constructor for NDTCore.
   * @author LegitYew0108(Wada Haruto)
   */
  explicit NDTCore();
  /**
   * @brief constructor with configs for NDTCore.
   * @author LegitYew0108(Wada Haruto)
   */
  explicit NDTCore(const NDTConfig &configs);
  /**
   * @brief destructor with configs for NDTCore.
   * @author LegitYew0108(Wada Haruto)
   */
  ~NDTCore();
  /**
   * @brief set all configurations for PureNDT3D
   * @author LegitYew0108(Wada Haruto)
   */
  void set_configurations(const NDTConfig &configs);
  /**
   * @brief add target points for PureNDT3D
   * @note If you want to replace all target points, use replace_target_points()
   * @author LegitYew0108(Wada Haruto)
   */
  void add_target_points(const std::vector<Point3D> &points);
  /**
   * @brief remove all target points for PureNDT3D
   * @author LegitYew0108(Wada Haruto)
   */
  void remove_target_points();
  /**
   * @brief add target points for PureNDT3D
   * @note If you want to add points to previously entered points, use
   * add_target_points().
   * @author LegitYew0108(Wada Haruto)
   */
  void replace_target_points(const std::vector<Point3D> &points);
  /**
   * @brief input points for PureNDT3D and returns aligned new_transform.
   * @author LegitYew0108(Wada Haruto)
   */
  Transform3D align(const std::vector<Point3D> &points,
                    const Transform3D &initial_transform);

protected:
  NDTConfig configs_;
  std::unique_ptr<VoxelGrid> voxel_grid_;
};

} // namespace PureNDT3D
