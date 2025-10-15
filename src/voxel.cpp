#include "voxel.hpp"
#include <unordered_map>

namespace PureNDT3D {
std::unordered_map<VoxelIndex, Voxel, VoxelHash> *VoxelGrid::get_voxels() {
  return &voxels_;
}
} // namespace PureNDT3D
