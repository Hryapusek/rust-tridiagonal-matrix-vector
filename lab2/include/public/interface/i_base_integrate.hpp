#pragma once

#include <defines.hpp>

#include <Eigen/Dense>

class IBaseIntegrate
{
 public:
  virtual auto integrate(
    Eigen::VectorXd const& start_v,
    Eigen::MatrixXd const& A,
    Eigen::VectorXd const& g,
    std::vector<Number_t> const& points
  ) -> Eigen::MatrixXd = 0;
};
