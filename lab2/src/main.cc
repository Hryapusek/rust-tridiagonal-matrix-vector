#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <defines.hpp>
#include <input_parameters.hpp>
#include <default_impl/main_matrix_calculator.hpp>
#include <interval_splitter.hpp>
#include <default_impl/euler_explicit_method.hpp>
#include <default_impl/euler_implicit_method.hpp>
#include <default_impl/rkf_method.hpp>

auto build_main_matrix(DefaultMainMatrixCalculator const& calc, double t)
  -> Eigen::SparseMatrix<Number_t>
{
  Eigen::SparseMatrix<Number_t> main_matrix(calc.r_points().size() - 1, calc.r_points().size() - 1);
  main_matrix.setZero();

  for(size_t row = 0; row < calc.r_points().size() - 1; ++row) {
    for(size_t col = 0; col < calc.r_points().size() - 1; ++col) {
      if(col == row) {
        main_matrix.insert(row, col) = calc.calc_c(row + 1, t);
      }
      else if(col == row + 1 and row != calc.r_points().size() - 1) {
        main_matrix.insert(row, col) = calc.calc_b(row + 1, t);
      }
      else if(col == row - 1 and row != 0) {
        main_matrix.insert(row, col) = calc.calc_a(row + 1, t);
      }
      else {
        (void)0;
      }
    }
  }

  return main_matrix;
}

auto build_g_vector(DefaultMainMatrixCalculator const& calc, Number_t t)
  -> Eigen::SparseVector<Number_t>
{
  auto g = Eigen::SparseVector<Number_t>(calc.r_points().size() - 1);
  for(size_t row = 0; row < calc.r_points().size() - 1; ++row) {
    g.insert(row) = calc.calc_g(row + 1, t);
  }
  return g;
}

auto integrate(
  DefaultMainMatrixCalculator const& calc,
  std::vector<Number_t> const& r_points,
  std::vector<Number_t> const& t_points,
  Eigen::VectorX<Number_t> const& start_v,
  IBaseIntegrate& method
) -> Eigen::MatrixX<Number_t>
{
  auto result = Eigen::MatrixXd(Eigen::MatrixXd::Zero(r_points.size() - 1, t_points.size()));
  result.col(0) = start_v;

  for(size_t i = 1; i < 2; ++i) {
    Eigen::SparseMatrix<Number_t> main_matrix;
    if (dynamic_cast<IEulerImplicitMethod*>(&method)) {
      main_matrix = build_main_matrix(calc, t_points.at(i));
    } else {
      main_matrix = build_main_matrix(calc, t_points.at(i - 1));
    }
    auto g = build_g_vector(calc, t_points.at(i));
    auto points = std::vector<Number_t>(t_points.cbegin() + i - 1, t_points.cbegin() + i + 1);
    auto integrated_result = method.integrate(result.col(i-1), main_matrix, g, points);
    result.col(i) = integrated_result.col(1);
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
  params->f = [](double r, double t) { return 3 * t + 6 * r + 1 - 2 / r; };

  auto r_points = split_interval(params->Rl, params->Rr, 1'000);
  DefaultMainMatrixCalculator calc(params, r_points);

  auto main_matrix = build_main_matrix(calc, 0);
  std::cout << main_matrix.topLeftCorner(5, 5) << std::endl;

  auto g_vector = build_g_vector(calc, 0);
  // std::cout << g_vector << std::endl;

  auto expected_func = [](double r, double t) { return t + 2 * r; };

  RKFMethod method;
  Eigen::SparseVector<Number_t> start_v(r_points.size() - 1);
  for(auto i = 0; i < r_points.size() - 1; ++i) {
    start_v.insert(i) = params->phi(r_points.at(i + 1));
  }
  auto t_points = split_interval(0, params->T, 10'000);
  auto result = integrate(calc, r_points, t_points, start_v, method);

  std::cout << "result: " << result.coeff(0, 1) << std::endl;
  std::cout << "expected: " << expected_func(r_points.at(1), t_points.at(1)) << std::endl;
}

int main()
{
  basic_example();
  return 0;
}
