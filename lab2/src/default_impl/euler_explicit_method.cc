#include <default_impl/euler_explicit_method.hpp>

#include <Eigen/Dense>

#include <contract/contract.hpp>

auto DefaultEulerExplicitMethod::integrate(
  Eigen::SparseVector<Number_t> const& start_v,
  Eigen::SparseMatrix<Number_t> const& A,
  Eigen::SparseVector<Number_t> const& g,
  std::vector<Number_t> const& points
) -> Eigen::SparseMatrix<Number_t>
{
  auto result = Eigen::SparseMatrix<Number_t>(
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

  result.col(0) = start_v;

  auto E = Eigen::SparseMatrix<Number_t>(A.rows(), A.rows());
  E.setIdentity();
  for(size_t i = 1; i < points.size(); ++i) {                 // Loop over intervals, not A.cols()
    auto H = points.at(i) - points.at(i - 1);                 // Step size
    result.col(i) = (E + H * A) * result.col(i - 1) + H * g;  // Euler update
  }

  return result;
}
