#pragma once

#include <interface/i_euler_explicit_method.hpp>

#include <Eigen/Dense>

class DefaultEulerExplicitMethod : public IEulerExplicitMethod
{
 public:
  auto integrate(
    Eigen::SparseVector<Number_t> const& start_v,
    Eigen::SparseMatrix<Number_t> const& A,
    Eigen::SparseVector<Number_t> const& g,
    std::vector<Number_t> const& points
  ) -> Eigen::SparseMatrix<Number_t> override;
};
