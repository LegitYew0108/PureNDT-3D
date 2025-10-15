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
using Transform3D = Eigen::Matrix3d;

enum class LogLevel { Error, Warning, Info, Debug };

using LoggerCallback = std::function<void(LogLevel, const std::string &)>;

class NDTCore {
public:
  /**
   * @brief default constructor for NDTCore.
   */
  explicit NDTCore();
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
  void remove_target_points();
  /**
   * @brief add target points for PureNDT3D
   * @note If you want to add points to previously entered points, use
   * add_target_points().
   */
  void replace_target_points(const std::vector<Point3D> &points);
  /**
   * @brief input points for PureNDT3D and returns aligned new_transform.
   */
  Transform3D align(const std::vector<Point3D> &points,
                    const Transform3D &initial_transform);

  /**
   * @brief setter for logger_
   */
  void set_logger(LoggerCallback logger);

protected:
  NDTConfig configs_;
  std::unique_ptr<VoxelGrid> voxel_grid_;

  void check_config(NDTConfig &configs_);

  LoggerCallback logger_;

  void log(LogLevel level, const std::string &message) {
    if (logger_) {
      logger_(level, message);
    }
  }
};

} // namespace PureNDT3D
