#include <iostream>
#include <Eigen/Dense>

#include <defines.hpp>
#include <input_parameters.hpp>
#include <default_impl/main_matrix_calculator.hpp>
#include <interval_splitter.hpp>
#include <default_impl/euler_explicit_method.hpp>
#include <default_impl/euler_implicit_method.hpp>

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

  auto interval = split_interval(params->Rl, params->Rr, 100);
  DefaultMainMatrixCalculator calc(params, interval);

  auto main_matrix = build_main_matrix(calc);
  // std::cout << main_matrix << std::endl;

  auto g_vector = build_g_vector(calc);
  // std::cout << g_vector << std::endl;

  DefaultEulerExplicitMethod method;
  Eigen::VectorXd start_v(interval.size() - 1);
  for(auto i = 0; i < interval.size() - 1; ++i) {
    start_v(i) = params->phi(interval.at(i + 1));
  }
  auto t_interval = split_interval(0, params->T, 100);
  auto result = method.integrate(start_v, main_matrix, g_vector, t_interval);

  std::cout << result(0, 1) << std::endl;
}

int main()
{
  basic_example();
  return 0;
}
