#include "ndt_core.hpp"
#include <gtest/gtest.h>

class NdtCoreFixture : public PureNDT3D::NDTCore {
public:
  PureNDT3D::NDTConfig get_config() { return configs_; }
};

TEST(NDTInterfaceTest, PositiveCase) {
  NdtCoreFixture ndt_core = NdtCoreFixture();

  PureNDT3D::NDTConfig config = PureNDT3D::NDTConfig();
  config.voxel_resolution_m_ = 0.1f;
  config.max_iterations = 2000;
  config.epsilon_trans = 0.001f;
  config.epsilon_rot = 0.1f;
  ndt_core.set_configurations(config);

  PureNDT3D::NDTConfig got_config = ndt_core.get_config();

  ASSERT_EQ(got_config.voxel_resolution_m_, 0.1f);
  ASSERT_EQ(got_config.max_iterations, 2000);
  ASSERT_EQ(got_config.epsilon_trans, 0.001f);
  ASSERT_EQ(got_config.epsilon_rot, 0.1f);
}
