#pragma once

#include <defines.hpp>

#include <Eigen/Dense>

class IEulerMethod
{
 public:
  virtual auto integrate(
    Eigen::VectorXd const& start_v,
    Eigen::MatrixXd& A,
    Eigen::VectorXd const& g,
    std::vector<Number_t> const& intervals
  ) -> Eigen::MatrixXd = 0;
};
