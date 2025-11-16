#pragma once

#include <functional>
#include <string>

namespace PureNDT3D {

enum class LogLevel { Error, Warning, Info, Debug };

using LoggerCallback = std::function<void(LogLevel, const std::string &)>;

struct NDTConfig {
  // === VoxelGridに関する設定===
  double voxel_resolution_m_{0.3};

  // === 最適化計算の設定===
  int max_iterations_{100};
  double score_threshold_{0.1};
  double levenberg_marquardt_lambda_{1e-6};

  double outlier_ratio_{0.05};
  double trans_alpha{0.1};
  double rot_alpha{0.05};

  bool use_second_order_derivative_;
  LoggerCallback logger_{nullptr};
};
} // namespace PureNDT3D
