#pragma once

#include <defines.hpp>

#include <Eigen/Sparse>

class IBaseIntegrate
{
 public:
  virtual auto integrate(
    Eigen::SparseVector<Number_t> const& start_v,
    Eigen::SparseMatrix<Number_t> const& A,
    Eigen::SparseVector<Number_t> const& g,
    std::vector<Number_t> const& points
  ) -> Eigen::SparseMatrix<Number_t> = 0;
};
