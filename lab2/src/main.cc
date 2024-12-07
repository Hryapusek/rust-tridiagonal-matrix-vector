#include <iostream>
#include <Eigen/Dense>

#include <defines.hpp>
#include <input_parameters.hpp>
#include <default_impl/main_matrix_calculator.hpp>
#include <interval_splitter.hpp>

/*
struct InputParameters {
  Number_t Rl;
  Number_t Rr; // [Rl, Rr]
  Number_t T; // [0, T]

  // First type condition
  T_Function_type v1; // u(rL) = v1(t)

  // Third type condition
  Number_t hi2; // -k * du/dr = hi2 * u(rR) - v2(t)
  R_Function_type phi;
  T_Function_type v2;

  // Just input functions
  R_T_Function_type k;
  R_T_Function_type q;
  R_T_Function_type f;
};
*/

auto build_main_matrix(DefaultMainMatrixCalculator const &calc) -> Eigen::MatrixXd {
  Eigen::MatrixXd main_matrix = Eigen::MatrixXd::Zero(calc.intervals().size() - 1, calc.intervals().size() - 1);
  for (size_t row = 0; row < calc.intervals().size() - 1; ++row) {
    for (size_t col = 0; col < calc.intervals().size() - 1; ++col) {
      if (col == row) {
        main_matrix(row, col) = calc.calc_c(row + 1);
      } else if (col == row + 1 and row != calc.intervals().size() - 1) {
        main_matrix(row, col) = calc.calc_b(row + 1);
      } else if (col == row - 1 and row != 0) {
        main_matrix(row, col) = calc.calc_a(row + 1);
      } else {
        (void)0;
      }
    }
  }

  return main_matrix;
}

void basic_example()
{
  std::shared_ptr<InputParameters> params = std::make_shared<InputParameters>();
  params->Rl = 1;
  params->Rr = 10;
  params->T = 10;
  params->v1 = [](double t) { return 2 + t; };
  params->hi2 = 3;
  params->phi = [](double r) { return 2 * r; };
  params->v2 = [](double t) { return 58 + 3 * t; };

  params->k = [](double r, double t) { return 1.0; };
  params->q = [](double r, double t) { return 3.0; };
  params->f = [](double r, double t) { return 6 * r + 3 * t - 2 / r + 1; };
  
  auto interval = split_interval(params->Rl, params->Rr, 10);
  DefaultMainMatrixCalculator calc(params, interval);

  auto main_matrix = build_main_matrix(calc);
  std::cout << main_matrix << std::endl;
}

int main()
{
  basic_example();
  return 0;
}
