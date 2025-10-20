#include "config.hpp"
#include "ndt_core.hpp"
#include <gtest/gtest.h>

class NdtCoreFixture : public PureNDT3D::NDTCore {
public:
  NdtCoreFixture() : PureNDT3D::NDTCore() {}
  NdtCoreFixture(PureNDT3D::NDTConfig &configs) : PureNDT3D::NDTCore(configs) {}
  PureNDT3D::NDTConfig get_config() { return configs_; }
};

TEST(NDTConfigTest, PositiveCase) {
  NdtCoreFixture ndt_core = NdtCoreFixture();

  PureNDT3D::NDTConfig config = PureNDT3D::NDTConfig();
  config.voxel_resolution_m_ = 0.1f;
  config.max_iterations_ = 2000;
  ndt_core.set_configurations(config);

  PureNDT3D::NDTConfig got_config = ndt_core.get_config();

  ASSERT_EQ(got_config.voxel_resolution_m_, 0.1f);
  ASSERT_EQ(got_config.max_iterations_, 2000);

  NdtCoreFixture ndt_core_configured = NdtCoreFixture(config);

  PureNDT3D::NDTConfig configured_config = ndt_core_configured.get_config();

  ASSERT_EQ(configured_config.voxel_resolution_m_, 0.1f);
  ASSERT_EQ(configured_config.max_iterations_, 2000);
}

TEST(NDTConfigTest, NegativeCase) {
  NdtCoreFixture ndt_core = NdtCoreFixture();
  PureNDT3D::NDTConfig config = PureNDT3D::NDTConfig();
  config.voxel_resolution_m_ = -0.1f;

  ASSERT_ANY_THROW(ndt_core.set_configurations(config));
}
