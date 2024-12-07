#include <default_impl/euler_explicit_method.hpp>

#include <gtest/gtest.h>

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <vector>
#include "default_impl/euler_implicit_method.hpp"

using namespace ::testing;

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <vector>
#include "default_impl/euler_explicit_method.hpp"
#include "default_impl/euler_implicit_method.hpp"  // Include the header for the implicit method

using namespace ::testing;

// Test for the implicit Euler integration method
TEST(EulerImplicitMethodTest, SimpleImplicitEulerIntegration) {
    // Set up inputs
    Eigen::MatrixXd A(2, 2);
    A << 0.1, 0.2,
         0.3, 0.4;

    Eigen::VectorXd g(2);
    g << 1.0, 2.0;

    Eigen::VectorXd start_v(2);
    start_v << 0.0, 0.0;

    std::vector<Number_t> intervals = {0.0, 0.1, 0.2, 0.3};  // Example intervals

    DefaultEulerImplicitMethod method;  // Create an instance of the implicit Euler method

    // Call the method
    Eigen::MatrixXd result = method.integrate(start_v, A, g, intervals);

    // Check the dimensions
    EXPECT_EQ(result.rows(), A.rows());
    EXPECT_EQ(result.cols(), intervals.size());

    // Expected results based on manual calculations (already precalculated for implicit Euler)
    Eigen::MatrixXd expected_result(2, 4);
    expected_result << 0.0, 0.095, 0.188, 0.2785,  // Implicit Euler results
                       0.0, 0.19, 0.374, 0.557;

    std::cout << result << std::endl;

    // Check if the result is approximately the expected result (we might adjust tolerance based on expected error)
    EXPECT_TRUE(result.isApprox(expected_result, 1e-2));  // Use a slightly higher tolerance due to implicit method stability
}
