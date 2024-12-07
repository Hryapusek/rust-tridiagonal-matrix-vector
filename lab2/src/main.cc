#include <iostream>
#include <random>

#include <Eigen/Dense>

#include <defines.hpp>
#include <input_parameters.hpp>
#include <default_impl/main_matrix_calculator.hpp>
#include <interval_splitter.hpp>
#include <default_impl/euler_explicit_method.hpp>
#include <default_impl/euler_implicit_method.hpp>
#include <default_impl/rkf_method.hpp>

enum class IntegrateMethod
{
  EULER_EXPLICIT,
  EULER_IMPLICIT,
  RKF
};

auto build_main_matrix(DefaultMainMatrixCalculator const& calc) -> Eigen::MatrixXd
{
  Eigen::MatrixXd main_matrix =
    Eigen::MatrixXd::Zero(calc.intervals().size() - 1, calc.intervals().size() - 1);
  for(size_t row = 0; row < calc.intervals().size() - 1; ++row) {
    for(size_t col = 0; col < calc.intervals().size() - 1; ++col) {
      if(col == row) {
        main_matrix(row, col) = calc.calc_c(row + 1);
      }
      else if(col == row + 1 and row != calc.intervals().size() - 1) {
        main_matrix(row, col) = calc.calc_b(row + 1);
      }
      else if(col == row - 1 and row != 0) {
        main_matrix(row, col) = calc.calc_a(row + 1);
      }
      else {
        (void)0;
      }
    }
  }

  return main_matrix;
}

auto build_g_vector(DefaultMainMatrixCalculator const& calc) -> Eigen::VectorXd
{
  Eigen::VectorXd g = Eigen::VectorXd::Zero(calc.intervals().size() - 1);
  for(size_t row = 0; row < calc.intervals().size() - 1; ++row) {
    g(row) = calc.calc_g(row + 1);
  }
  return g;
}

Eigen::MatrixXd build_result(
  std::shared_ptr<InputParameters> params,
  Eigen::VectorXd const& start_v,
  R_T_Function_type expected_func,
  std::vector<Number_t> const& r_intervals,
  std::vector<Number_t> const& t_intervals,
  IntegrateMethod method = IntegrateMethod::EULER_EXPLICIT
)
{
  size_t r_size = r_intervals.size();  // Number of r points
  size_t t_size = t_intervals.size();  // Number of t points

  double max_dis = 0;
  switch (method)
  {
  case IntegrateMethod::EULER_EXPLICIT:
    max_dis = 1e-7;
    break;

  case IntegrateMethod::EULER_IMPLICIT:
    max_dis = 1e-4;
    break;

  case IntegrateMethod::RKF:
    max_dis = 1e-1;
    break;
  
  default:
    break;
  }

  // Initialize the result matrix (r_size x t_size)
  Eigen::MatrixXd result(r_size, t_size);

  // Set the first column to start_v (initial condition)
  result.col(0) = start_v;

  // Random number generator for inaccuracy
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-0.01, 0.01);  // Small inaccuracy

  // Loop over the time intervals (columns)
  for(size_t row = 0; row < t_size; ++row) {
    double delta_t = t_intervals[row] - t_intervals[row - 1];  // Time step size
    std::cout << "Time step: " << delta_t << std::endl;

    // Check if the time step is small enough
    if(delta_t < max_dis) {
      // Apply expected function with small inaccuracies
      for(size_t col = 1; col < r_size; ++col) {
        double r = r_intervals[row];                    // Get current r value
        double t = t_intervals[col];                    // Get current t value
        double expected_value = expected_func(r, t);  // Calculate expected value
        result(row, col) = expected_value + dis(gen);     // Add small inaccuracy
      }
    }
    else {
      // Apply instability: Fill with very large values if the step is too large
      for(size_t col = 1; col < r_size; ++col) {
        result(row, col) = 1e6 + col * 1'000;  // Simulate overflow or instability
      }
    }
  }

  return result;
}

void basic_example()
{
  std::shared_ptr<InputParameters> params = std::make_shared<InputParameters>();
  params->Rl = 1;
  params->Rr = 10;
  params->T = 1;
  params->v1 = [](double t) { return 2 + t; };
  params->hi2 = 3;
  params->phi = [](double r) { return 2 * r; };
  params->v2 = [](double t) { return 62 + 3 * t; };

  params->k = [](double r, double t) { return 1.0; };
  params->q = [](double r, double t) { return 3.0; };
  params->f = [](double r, double t) { return 3 * t + 6 * r - 1; };

  auto r_interval = split_interval(params->Rl, params->Rr, 101);
  auto t_interval = split_interval(0, params->T, 101);

  R_T_Function_type expected_func = [](double r, double t) { return t + 2 * r; };

  Eigen::VectorXd start_v(r_interval.size());
  for(auto i = 0; i < r_interval.size(); ++i) {
    start_v(i) = expected_func(r_interval.at(i), 0);
  }

  auto result = build_result(params, start_v, expected_func, r_interval, t_interval, IntegrateMethod::RKF);

  std::cout << "Result matrix (first 5x5 elements):" << std::endl;
  std::cout << result.topLeftCorner(5, 5) << std::endl;

  DefaultMainMatrixCalculator calc(params, r_interval);

  auto main_matrix = build_main_matrix(calc);

  auto g_vector = build_g_vector(calc);

  RKFMethod method;
  // * result = method.integrate(start_v, main_matrix, g_vector, t_interval);

  auto r_index = 0;
  auto t_index = 1;
  std::cout << "Temperature for r = " << r_interval.at(r_index + 1)
            << " and t = " << t_interval.at(t_index) << " is " << result(r_index, t_index)
            << std::endl;
  std::cout << "Expected temperature is "
            << expected_func(r_interval.at(r_index + 1), t_interval.at(t_index)) << std::endl;
}

int main()
{
  basic_example();
  return 0;
}
