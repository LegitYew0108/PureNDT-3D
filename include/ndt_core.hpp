#pragma once

#include "config.hpp"
#include "eigen3/Eigen/Core"
#include "ndt_matcher.hpp"
#include "voxel.hpp"
#include <functional>
#include <memory>
#include <vector>

namespace PureNDT3D {

using Point3D = Eigen::Vector3d;
using Covariance3D = Eigen::Matrix3d;
using Transform4D = Eigen::Matrix4d;

void log(LoggerCallback logger, LogLevel level, const char *message, ...);

class NDTCore {
public:
  /**
   * @brief constructor with configs for NDTCore.
   */
  explicit NDTCore(const NDTConfig &configs);
  /**
   * @brief destructor with configs for NDTCore.
   */
  ~NDTCore();
  /**
   * @brief set all configurations for PureNDT3D
   */
  void set_configurations(const NDTConfig &configs);
  /**
   * @brief add target points for PureNDT3D
   * @note If you want to replace all target points, use replace_target_points()
   */
  void add_target_points(const std::vector<Point3D> &points);
  /**
   * @brief remove all target points for PureNDT3D
   */
  void remove_all_target_points();
  /**
   * @brief add target points for PureNDT3D
   * @note If you want to add points to previously entered points, use
   * add_target_points().
   */
  void replace_target_points(const std::vector<Point3D> &points);
  /**
   * @brief input points for PureNDT3D and returns aligned new_transform.
   */
  Transform4D align(const std::vector<Point3D> &points,
                    const Transform4D &initial_transform);

protected:
  NDTConfig configs_;
  std::unique_ptr<VoxelGrid> voxel_grid_;
  std::unique_ptr<NDTMatcher> matcher_;

  void check_config(NDTConfig &configs_);
};

} // namespace PureNDT3D
