#include <default_impl/rkf_method.hpp>

#include <Eigen/Sparse>

auto RKFMethod::integrate(
  Eigen::VectorX<Number_t> const& start_v,
  Eigen::MatrixX<Number_t> const& A,
  Eigen::SparseVector<Number_t> const& g,
  std::vector<Number_t> const& points
) -> Eigen::SparseMatrix<Number_t>
{
  auto result = Eigen::MatrixX<Number_t>(A.rows(), A.cols());
  result.col(0) = start_v.sparseView();

  for(size_t i = 1; i < points.size(); ++i) {
    auto H = points.at(i) - points.at(i - 1);

    Eigen::MatrixX<Number_t> u = result.col(i - 1);

    Eigen::MatrixX<Number_t> k1 = (H * (A * u + g));
    Eigen::MatrixX<Number_t> k2 = (H * (A * (u + 0.5 * k1) + g));
    Eigen::MatrixX<Number_t> k3 = (H * (A * (u + 0.5 * k2) + g));
    Eigen::MatrixX<Number_t> k4 = (H * (A * (u + k3) + g));

    result.col(i) = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  }

  return result.sparseView();
}
