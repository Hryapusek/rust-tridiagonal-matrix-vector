#pragma once

#include <interface/i_euler_implicit_method.hpp>

#include <Eigen/Dense>

class DefaultEulerImplicitMethod : public IEulerImplicitMethod
{
 public:
  auto integrate(
    Eigen::VectorXd const& start_v,
    Eigen::MatrixXd& A,
    Eigen::VectorXd const& g,
    std::vector<Number_t> const& intervals
  ) -> Eigen::MatrixXd override;
};
