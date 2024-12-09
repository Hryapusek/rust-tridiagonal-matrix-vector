#pragma once

#include <defines.hpp>

#include <Eigen/Sparse>

class IBaseIntegrate
{
 public:
  virtual auto integrate(
    Eigen::VectorX<Number_t> const& start_v,
    Eigen::MatrixX<Number_t> const& A,
    Eigen::SparseVector<Number_t> const& g,
    std::vector<Number_t> const& points
  ) -> Eigen::SparseMatrix<Number_t> = 0;
};
