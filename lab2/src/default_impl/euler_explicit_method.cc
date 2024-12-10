#include <default_impl/euler_explicit_method.hpp>

#include <Eigen/Dense>

#include <contract/contract.hpp>

auto DefaultEulerExplicitMethod::integrate(
  Eigen::VectorX<Number_t> const& start_v,
  Eigen::MatrixX<Number_t> const& A,
  Eigen::SparseVector<Number_t> const& g,
  std::vector<Number_t> const& points
) -> Eigen::SparseMatrix<Number_t>
{
  auto result = Eigen::MatrixX<Number_t>(
    A.rows(),
    points.size()
  );  // Adjust columns based on intervals size

  contract(fun)
  {
    precondition(A.rows() == A.cols(), "A must be a square matrix");
    precondition(A.rows() == start_v.rows(), "A and start_v must have the same number of rows");
    precondition(A.rows() == g.rows(), "A and g must have the same number of rows");
    postcondition(
      result.cols() == points.size(),
      "result must have the same number of columns as intervals"
    );
  };

  result.col(0) = start_v.sparseView();

  auto E = Eigen::SparseMatrix<Number_t>(A.rows(), A.rows());
  E.setIdentity();
  auto H = points.at(1) - points.at(0);                 // Step size
  for(size_t i = 1; i < 2; ++i) {                 // Loop over intervals, not A.cols()
  // for(size_t i = 1; i < points.size(); ++i) {
    result.col(i) = (E + H * A) * result.col(i - 1) + H * g;  // Euler update
  }

  return result.sparseView();
}
