#include "ndt_core.hpp"
#include <gtest/gtest.h>

TEST(NDTInterfaceTest, PositiveCase) {
  PureNDT3D::NDTCore ndt_core = PureNDT3D::NDTCore();
  PureNDT3D::NDTConfig config = PureNDT3D::NDTConfig();
  ndt_core.set_configurations(config);
}
