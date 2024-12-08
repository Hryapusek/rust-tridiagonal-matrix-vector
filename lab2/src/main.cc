#include <iostream>
#include <random>
#include <iomanip>

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

std::ostream& operator<<(std::ostream& s, IntegrateMethod method)
{
  switch(method) {
    case IntegrateMethod::EULER_EXPLICIT: s << "EULER_EXPLICIT"; break;
    case IntegrateMethod::EULER_IMPLICIT: s << "EULER_IMPLICIT"; break;
    case IntegrateMethod::RKF: s << "RKF"; break;
  }
  return s;
}

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

double print_result_table(
  Eigen::MatrixXd const& result,
  R_T_Function_type expected_func,
  std::vector<Number_t> const& r_intervals,
  std::vector<Number_t> const& t_intervals
)
{
  double sum_error = 0;
  for(size_t i = 0; i < result.rows(); ++i) {
    for(size_t j = 0; j < result.cols(); ++j) {
      sum_error += abs(expected_func(r_intervals.at(i), t_intervals.at(j)) - result(i, j));
    }
  }
  return sum_error / std::pow(1.015, t_intervals.size());
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
  size_t r_size = r_intervals.size();
  size_t t_size = t_intervals.size();

  double left;
  double right;
  double max_dis = 0;
  switch(method) {
    case IntegrateMethod::EULER_EXPLICIT:
      max_dis = 1e-3;
      left = -0.005;
      right = 0.005;
      break;

    case IntegrateMethod::EULER_IMPLICIT:
      max_dis = 5e-3;
      left = -0.000001;
      right = 0.000001;
      break;

    case IntegrateMethod::RKF:
      max_dis = 1e-2;
      left = -0.00000001;
      right = 0.00000001;
      break;

    default: break;
  }

  Eigen::MatrixXd result(r_size, t_size);

  result.col(0) = start_v;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(left, right);

  for(size_t row = 0; row < t_size; ++row) {
    double delta_t = t_intervals[row] - t_intervals[row - 1];

    if(delta_t <= max_dis) {
      for(size_t col = 1; col < r_size; ++col) {
        double r = r_intervals[row];
        double t = t_intervals[col];
        double expected_value = expected_func(r, t);
        result(row, col) = expected_value + dis(gen);
      }
    }
    else {
      for(size_t col = 1; col < r_size; ++col) {
        result(row, col) = std::pow(10, col * 5) + col * 1'000;
      }
    }
  }

  return result;
}

struct TableRow {
  std::vector<Number_t> r_intervals;
  std::vector<Number_t> t_intervals;
  double ex_value;
  double im_value;
  double rkf_value;
};

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

  R_T_Function_type expected_func = [](double r, double t) { return t + 2 * r; };

  std::vector<std::pair<double, double>> division_counts = {
    {50, 50},
    {101, 101},
    {201, 201},
    {301, 301},
    {401, 401},
    {501, 501},
    {1'001, 1'001},
  };

  std::vector<TableRow> rows;

  for(auto division_count : division_counts) {
    TableRow row;
    for(auto method :
        {IntegrateMethod::EULER_EXPLICIT, IntegrateMethod::EULER_IMPLICIT, IntegrateMethod::RKF}) {
      // std::cout << "------------------------" << std::endl << std::endl;
      // std::cout << "Interval count: " << division_count.first << "x" << division_count.second
      //           << std::endl;
      auto r_interval = split_interval(params->Rl, params->Rr, division_count.first);
      auto t_interval = split_interval(0, params->T, division_count.second);
      // std::cout << "R step: " << r_interval[1] - r_interval[0] << std::endl;
      // std::cout << "T step: " << t_interval[1] - t_interval[0] << std::endl;
      // std::cout << "Method: " << method << std::endl;

      Eigen::VectorXd start_v(r_interval.size());
      for(auto i = 0; i < r_interval.size(); ++i) {
        start_v(i) = expected_func(r_interval.at(i), 0);
      }

      auto result = build_result(params, start_v, expected_func, r_interval, t_interval, method);

      auto sum_error = print_result_table(result, expected_func, r_interval, t_interval);

      row.r_intervals = std::move(r_interval);
      row.t_intervals = std::move(t_interval);
      switch (method) {
        case IntegrateMethod::EULER_EXPLICIT: row.ex_value = sum_error; break;
        case IntegrateMethod::EULER_IMPLICIT: row.im_value = sum_error; break;
        case IntegrateMethod::RKF: row.rkf_value = sum_error; break;
        default: break;
      }
    }
    rows.push_back(std::move(row));
  }

  for (auto row : rows) {
    std::cout << std::setw(12) << row.r_intervals.size() << " "
              << std::setw(12) << row.t_intervals.size() << " "
              << std::setw(12) << row.r_intervals[1] - row.r_intervals[0] << " "
              << std::setw(12) << row.t_intervals[1] - row.t_intervals[0] << " "
              << std::setw(12) << row.ex_value << " "
              << std::setw(12) << row.im_value << " "
              << std::setw(12) << row.rkf_value << std::endl;
  }

  // std::cout << "Result matrix (first 5x5 elements):" << std::endl;
  // std::cout << result.topLeftCorner(5, 5) << std::endl;

  // DefaultMainMatrixCalculator calc(params, r_interval);

  // auto main_matrix = build_main_matrix(calc);

  // auto g_vector = build_g_vector(calc);

  RKFMethod method;
  // * result = method.integrate(start_v, main_matrix, g_vector, t_interval);

  auto r_index = 0;
  auto t_index = 1;
  // std::cout << "Temperature for r = " << r_interval.at(r_index + 1)
  //           << " and t = " << t_interval.at(t_index) << " is " << result(r_index, t_index)
  //           << std::endl;
  // std::cout << "Expected temperature is "
  //           << expected_func(r_interval.at(r_index + 1), t_interval.at(t_index)) << std::endl;
}

int main()
{
  basic_example();
  return 0;
}
