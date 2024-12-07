#include <default_impl/euler_explicit_method.hpp>

#include <Eigen/Dense>

#include <contract/contract.hpp>

auto DefaultEulerExplicitMethod::integrate(
  Eigen::VectorXd const& start_v,
  Eigen::MatrixXd& A,
  Eigen::VectorXd const& g,
  std::vector<Number_t> const& intervals
) -> Eigen::MatrixXd
{
  Eigen::MatrixXd result =
    Eigen::MatrixXd::Zero(A.rows(), intervals.size());  // Adjust columns based on intervals size

  contract(fun)
  {
    precondition(A.rows() == A.cols(), "A must be a square matrix");
    precondition(A.rows() == start_v.rows(), "A and start_v must have the same number of rows");
    precondition(A.rows() == g.rows(), "A and g must have the same number of rows");
    postcondition(
      result.cols() == intervals.size(),
      "result must have the same number of columns as intervals"
    );
  };

  result.col(0) = start_v;

  auto E = Eigen::MatrixXd::Identity(A.rows(), A.rows());
  for(size_t i = 1; i < intervals.size(); ++i) {              // Loop over intervals, not A.cols()
    auto H = intervals.at(i) - intervals.at(i - 1);           // Step size
    result.col(i) = (E + H * A) * result.col(i - 1) + H * g;  // Euler update
  }

  return result;
}
