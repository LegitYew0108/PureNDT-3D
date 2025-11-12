#include "eigen3/Eigen/Core"

using Transform4D = Eigen::Matrix4d;
using TransformVec6D = Eigen::Vector<double, 6>;

struct TransformDatas {
  Transform4D mat;
  TransformVec6D vec;
};
