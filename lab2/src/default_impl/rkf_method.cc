#include <default_impl/rkf_method.hpp>

auto RKFMethod::integrate(
  Eigen::VectorXd const& start_v,
  Eigen::MatrixXd& A,
  Eigen::VectorXd const& g,
  std::vector<Number_t> const& intervals
) -> Eigen::MatrixXd
{
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(A.rows(), A.cols());
  result.col(0) = start_v;
  
  for (size_t i = 1; i < A.cols() - 1; ++i) {
    auto H = intervals.at(i) - intervals.at(i - 1);
    
    Eigen::VectorXd u = result.col(i - 1);
    
    Eigen::VectorXd k1 = H * (A * u + g);
    Eigen::VectorXd k2 = H * (A * (u + 0.5 * k1) + g);
    Eigen::VectorXd k3 = H * (A * (u + 0.5 * k2) + g);
    Eigen::VectorXd k4 = H * (A * (u + k3) + g);
    
    result.col(i) = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
  }
  
  return result;
}
