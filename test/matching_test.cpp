#include "config.hpp"
#include "ndt_core.hpp"
#include "transform.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <memory>

class NdtMatchingTest : public testing::Test {
protected:
  std::unique_ptr<PureNDT3D::NDTCore> core_;
  virtual void SetUp() {
    PureNDT3D::NDTConfig config;
    config.voxel_resolution_m_ = 0.5;
    config.logger_ = setupLogger();
    config.max_iterations_ = 20;

    core_ = std::make_unique<PureNDT3D::NDTCore>(config);
  }
  PureNDT3D::LoggerCallback setupLogger() {
    return [this](PureNDT3D::LogLevel level, const std::string &message) {
      switch (level) {
      case PureNDT3D::LogLevel::Debug:
        std::cout << "[Debug][PureNDT3Dlib] " << message << std::endl;
        break;
      case PureNDT3D::LogLevel::Info:
        std::cout << "[Info][PureNDT3Dlib] " << message << std::endl;
        break;
      case PureNDT3D::LogLevel::Warning:
        std::cout << "[Warn][PureNDT3Dlib] " << message << std::endl;
        break;
      case PureNDT3D::LogLevel::Error:
        std::cout << "[Error][PureNDT3Dlib] " << message << std::endl;
        break;
      }
    };
  }
};

TEST_F(NdtMatchingTest, OneVoxel) {
  std::vector<Point3D> points_vec;
  points_vec.emplace_back(0.5, 0.5, 0.5);
  points_vec.emplace_back(0.6, 0.6, 0.5);
  points_vec.emplace_back(0.5, 0.8, 0.5);
  points_vec.emplace_back(0.7, 0.9, 0.5);
  points_vec.emplace_back(0.3, 0.5, 0.5);
  points_vec.emplace_back(0.5, 0.5, 0.2);
  points_vec.emplace_back(0.6, 0.6, 0.2);
  points_vec.emplace_back(0.5, 0.8, 0.2);
  points_vec.emplace_back(0.7, 0.9, 0.2);
  points_vec.emplace_back(0.3, 0.5, 0.2);
  points_vec.emplace_back(0.5, 0.5, 0.7);
  points_vec.emplace_back(0.6, 0.6, 0.7);
  points_vec.emplace_back(0.5, 0.8, 0.7);
  points_vec.emplace_back(0.7, 0.9, 0.7);
  points_vec.emplace_back(0.3, 0.5, 0.7);
  core_->add_target_points(points_vec);
  TransformVec6D t_vec = TransformVec6D::Zero();
  PureNDT3D::TransformType t(t_vec);
  PureNDT3D::TransformType zero_transform = core_->align(points_vec, t);
  TransformVec6D exp_zero_vec = zero_transform.get_vector();
  std::cout << "exp zero vector: (" << exp_zero_vec(0) << ", "
            << exp_zero_vec(1) << ", " << exp_zero_vec(2) << ", "
            << exp_zero_vec(3) << ", " << exp_zero_vec(4) << ", "
            << exp_zero_vec(5) << ")" << std::endl;
  EXPECT_EQ(exp_zero_vec(0), 0.0);
  EXPECT_EQ(exp_zero_vec(1), 0.0);
  points_vec.clear();
  points_vec.emplace_back(0.3, 0.3, 0.5);
  points_vec.emplace_back(0.4, 0.4, 0.5);
  points_vec.emplace_back(0.3, 0.6, 0.5);
  points_vec.emplace_back(0.5, 0.7, 0.5);
  points_vec.emplace_back(0.3, 0.3, 0.5);
  points_vec.emplace_back(0.3, 0.3, 0.2);
  points_vec.emplace_back(0.4, 0.4, 0.2);
  points_vec.emplace_back(0.3, 0.6, 0.2);
  points_vec.emplace_back(0.5, 0.7, 0.2);
  points_vec.emplace_back(0.1, 0.3, 0.2);
  points_vec.emplace_back(0.3, 0.3, 0.7);
  points_vec.emplace_back(0.4, 0.4, 0.7);
  points_vec.emplace_back(0.3, 0.6, 0.7);
  points_vec.emplace_back(0.5, 0.7, 0.7);
  points_vec.emplace_back(0.1, 0.3, 0.7);
  PureNDT3D::TransformType new_transform = core_->align(points_vec, t);
  TransformVec6D result_vec = new_transform.get_vector();
  std::cout << "result was:" << result_vec(0) << result_vec(1) << result_vec(2)
            << result_vec(3) << result_vec(4) << result_vec(5) << std::endl;
  ASSERT_GT(result_vec(0), 0.0);
  ASSERT_GT(result_vec(1), 0.0);
}
