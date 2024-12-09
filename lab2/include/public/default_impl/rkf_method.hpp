#pragma once

#include <interface/i_base_integrate.hpp>

#include <Eigen/Dense>

class RKFMethod : public IBaseIntegrate
{
 public:
  auto integrate(
    Eigen::VectorX<Number_t> const& start_v,
    Eigen::MatrixX<Number_t> const& A,
    Eigen::SparseVector<Number_t> const& g,
    std::vector<Number_t> const& points
  ) -> Eigen::SparseMatrix<Number_t> override;
};
