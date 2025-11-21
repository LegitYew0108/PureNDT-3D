#include "config.hpp"
#include "ndt_core.hpp"
#include "transform.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <vector>

class NdtMatchingTest : public testing::Test {
protected:
  std::unique_ptr<PureNDT3D::NDTCore> core_;

  void SetUp() override {
    PureNDT3D::NDTConfig config;
    config.voxel_resolution_m_ = 1.0;
    config.logger_ = setupLogger();
    config.max_iterations_ = 20;
    config.trans_alpha = 0.01;
    config.rot_alpha = 0.005;

    core_ = std::make_unique<PureNDT3D::NDTCore>(config);
  }

  PureNDT3D::LoggerCallback setupLogger() {
    return [this](PureNDT3D::LogLevel level, const std::string &message) {
      if (level == PureNDT3D::LogLevel::Info) {
        std::cout << "[PureNDT3Dlib][Info] " << message << "\n" << std::endl;
      }
      if (level == PureNDT3D::LogLevel::Debug) {
        std::cout << "[PureNDT3Dlib][Debug] " << message << "\n" << std::endl;
      }
      if (level == PureNDT3D::LogLevel::Warning) {
        std::cout << "[PureNDT3Dlib][Warning] " << message << "\n" << std::endl;
      }
      if (level == PureNDT3D::LogLevel::Error) {
        std::cout << "[PureNDT3Dlib][Error] " << message << "\n" << std::endl;
      }
    };
  }

  // Helper to generate a grid of points
  std::vector<PureNDT3D::Point3D>
  createGridPoints(double start_x, double start_y, double start_z, int count_x,
                   int count_y, int count_z, double step) {
    std::vector<PureNDT3D::Point3D> points;
    points.reserve(count_x * count_y * count_z);
    for (int i = 0; i < count_x; ++i) {
      for (int j = 0; j < count_y; ++j) {
        for (int k = 0; k < count_z; ++k) {
          points.emplace_back(start_x + i * step, start_y + j * step,
                              start_z + k * step);
        }
      }
    }
    return points;
  }

  void printPoints(std::vector<PureNDT3D::Point3D> points) {
    for (const auto &point : points) {
      std::cout << "Point: (" << point.x() << ", " << point.y() << ", "
                << point.z() << ")" << std::endl;
    }
  }

  void printTransform(PureNDT3D::TransformType transform) {
    PureNDT3D::TransformVec6D vec = transform.get_vector();
    std::cout << "Transform: (" << vec(0) << " ," << vec(1) << " ," << vec(2)
              << " ," << vec(3) << " ," << vec(4) << " ," << vec(5) << ")"
              << std::endl;
  }
};

TEST_F(NdtMatchingTest, ZeroTransform) {
  // Create a target point cloud
  auto target_points = createGridPoints(0.3, 0.5, 0.2, 3, 3, 3, 0.1);
  core_->add_target_points(target_points);

  // Initial transform is zero
  Eigen::Vector<double, 6> t_vec = Eigen::Vector<double, 6>::Zero();
  PureNDT3D::TransformType initial_transform(t_vec);

  // Align the same points (source = target)
  PureNDT3D::TransformType result_transform =
      core_->align(target_points, initial_transform);
  printTransform(result_transform);
  PureNDT3D::TransformVec6D result_vec = result_transform.get_vector();

  // Expect the result to be very close to zero
  for (int i = 0; i < 6; ++i) {
    EXPECT_NEAR(result_vec(i), 0.0, 1e-6)
        << "Element " << i << " should be near 0";
  }
}

TEST_F(NdtMatchingTest, ShiftedTransform) {
  // Create a target point cloud
  auto target_points = createGridPoints(0.3, 0.5, 0.2, 4, 1, 4, 0.2);
  core_->add_target_points(target_points);

  // Create a source point cloud shifted by (-0.2, 0.0, 0.0)
  // So the true transform from source to target should be (+0.2, +0.0, +0.0)
  auto source_points = createGridPoints(0.1, 0.5, 0.2, 4, 1, 4, 0.2);

  TransformVec6D t_vec = TransformVec6D::Zero();
  PureNDT3D::TransformType initial_transform(t_vec);

  PureNDT3D::TransformType result_transform =
      core_->align(source_points, initial_transform);
  TransformVec6D result_vec = result_transform.get_vector();
  printTransform(result_transform);

  // We expect translation x and y to be around 0.2
  EXPECT_NEAR(result_vec(0), 0.2, 0.05); // x translation
  EXPECT_NEAR(result_vec(1), 0.0, 0.05); // y translation
  EXPECT_NEAR(result_vec(2), 0.0, 0.05); // z translation

  // Rotations should be near 0
  EXPECT_NEAR(result_vec(3), 0.0, 0.05); // roll
  EXPECT_NEAR(result_vec(4), 0.0, 0.05); // pitch
  EXPECT_NEAR(result_vec(5), 0.0, 0.05); // yaw
}
