#pragma once

namespace PureNDT3D {
struct NDTConfig {
  // === VoxelGridに関する設定===
  double voxel_resolution_m_{0.3f};

  // === 最適化計算の設定===
  int max_iterations{30};
  double epsilon_trans{1e-2};
  double epsilon_rot{1e-2};
};
} // namespace PureNDT3D
