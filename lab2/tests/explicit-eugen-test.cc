#include <default_impl/euler_explicit_method.hpp>

#include <gtest/gtest.h>

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <vector>
#include "default_impl/euler_explicit_method.hpp"  // Include your header for DefaultEulerExplicitMethod

using namespace ::testing;

// Test for a standard integration case
TEST(EulerExplicitMethodTest, StandardIntegration)
{
  // Set up inputs
  Eigen::MatrixXd A(2, 2);
  A << 0.1, 0.2, 0.3, 0.4;

  Eigen::VectorXd g(2);
  g << 1.0, 2.0;

  Eigen::VectorXd start_v(2);
  start_v << 0.0, 0.0;

  std::vector<Number_t> intervals = {0.0, 0.1, 0.2, 0.3};  // Example intervals

  DefaultEulerExplicitMethod method;

  // Call the method
  Eigen::MatrixXd result = method.integrate(start_v, A, g, intervals);

  // Check the result dimensions
  EXPECT_EQ(result.rows(), A.rows());
  EXPECT_EQ(result.cols(), intervals.size());

  // Check the first column (which should equal start_v)
  EXPECT_TRUE(result.col(0).isApprox(start_v));

  // You can also check further steps based on your expected output
  // Here you should add checks based on known outcomes for your method
}

// Test for edge case: single interval
TEST(EulerExplicitMethodTest, SingleInterval)
{
  // Set up inputs for a single interval (this is an edge case)
  Eigen::MatrixXd A(2, 2);
  A << 0.1, 0.2, 0.3, 0.4;

  Eigen::VectorXd g(2);
  g << 1.0, 2.0;

  Eigen::VectorXd start_v(2);
  start_v << 0.0, 0.0;

  std::vector<Number_t> intervals = {0.0, 0.1};  // Only one step

  DefaultEulerExplicitMethod method;

  // Call the method
  Eigen::MatrixXd result = method.integrate(start_v, A, g, intervals);

  // Check dimensions
  EXPECT_EQ(result.rows(), A.rows());
  EXPECT_EQ(result.cols(), intervals.size());

  // Check that the result for the first step is correctly calculated
  EXPECT_TRUE(result.col(0).isApprox(start_v));
}

TEST(EulerExplicitMethodTest, SimpleEulerIntegration)
{
  // Set up inputs
  Eigen::MatrixXd A(2, 2);
  A << 0.1, 0.2, 0.3, 0.4;

  Eigen::VectorXd g(2);
  g << 1.0, 2.0;

  Eigen::VectorXd start_v(2);
  start_v << 0.0, 0.0;

  std::vector<Number_t> intervals = {0.0, 0.1, 0.2, 0.3};  // Example intervals

  DefaultEulerExplicitMethod method;

  // Call the method
  Eigen::MatrixXd result = method.integrate(start_v, A, g, intervals);

  // Check the dimensions
  EXPECT_EQ(result.rows(), A.rows());
  EXPECT_EQ(result.cols(), intervals.size());

  // Expected results based on manual calculations
  Eigen::MatrixXd expected_result(2, 4);
  expected_result << 0.0, 0.1, 0.205, 0.316, 0.0, 0.2, 0.411, 0.626;

  std::cout << result << std::endl;

  // Check if the result is approximately the expected result
  EXPECT_TRUE(result.isApprox(expected_result, 1e-4)) << "Result does not match expected result";
}
