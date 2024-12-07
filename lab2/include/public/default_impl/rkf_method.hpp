#pragma once

#include <interface/i_euler_explicit_method.hpp>

#include <Eigen/Dense>

class RKFMethod
{
 public:
  auto integrate(
    Eigen::VectorXd const& start_v,
    Eigen::MatrixXd& A,
    Eigen::VectorXd const& g,
    std::vector<Number_t> const& intervals
  ) -> Eigen::MatrixXd;
};
