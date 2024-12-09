#include <default_impl/main_matrix_calculator.hpp>

#include <cassert>

#include <contract/contract.hpp>

#include <interval_splitter.hpp>

auto DefaultMainMatrixCalculator::calc_a(size_t index) const -> Number_t
{
  contract(fun)
  {
    precondition(
      index != 0,
      "You should not calculate anything for index == 0 - you already have v function"
    );
    precondition(index != 1, "You should not calculate a for index == 1");
    precondition(index < r_points_.size(), "index out of range");
  };
  auto mid_point = middle_point(r_points_, index);
  auto k_value = params_->k(mid_point, 0);
  auto h_value = calc_h(r_points_, index);

  return mid_point * k_value / h_value;
}

auto DefaultMainMatrixCalculator::calc_b(size_t index) const -> Number_t
{
  contract(fun)
  {
    precondition(
      index != 0,
      "You should not calculate anything for index == 0 - you already have v function"
    );
    precondition(
      index != r_points_.size() - 1,
      "You should not calculate b for index == intervals_.size() - 1"
    );
    precondition(index < r_points_.size(), "index out of range");
  };
  auto up = middle_point(r_points_, index + 1) * params_->k(middle_point(r_points_, index + 1), 0);
  auto down = calc_h(r_points_, index + 1);
  return up / down;
}

/* clang-format off */
auto DefaultMainMatrixCalculator::calc_c(size_t index) const -> Number_t
{
  contract(fun) {
    precondition(index != 0, "You should not calculate anything for index == 0 - you already have v function");
    precondition(index < r_points_.size(), "index out of range");
  };
  if (index == r_points_.size() - 1) {
    return - middle_point(r_points_, index) * params_->k(middle_point(r_points_, index), 0)
                 / calc_h(r_points_, index)
            - params_->hi2
            - params_->q(r_points_[index], 0)
                 * calc_cross_h(r_points_, index);
  } else {
    auto m_point_index = middle_point(r_points_, index);
    auto m_point_index_plus_1 = middle_point(r_points_, index + 1);
    auto k_index = params_->k(m_point_index, 0);
    auto k_index_plus_1 = params_->k(m_point_index_plus_1, 0);
    auto h_index_plus_1 = calc_h(r_points_, index + 1);
    auto cross_h_index = calc_cross_h(r_points_, index);
    auto r_index = r_points_[index];
    auto q_index = params_->q(r_index, 0);

    return -m_point_index * k_index / h_index_plus_1
           -m_point_index_plus_1 * k_index_plus_1 / h_index_plus_1
           -q_index * calc_cross_h(r_points_, index);
  }
  assert(false);
}

/* clang-format on */

/* clang-format off */
auto DefaultMainMatrixCalculator::calc_g(size_t index) const -> Number_t {
  contract(fun) {
    precondition(index != 0, "You should not calculate anything for index == 0 - you already have v function");
  };
  if (index == 1) {
    return params_->f(r_points_[index], 0) +
               middle_point(r_points_, index) 
               * params_->k(middle_point(r_points_, index), 0)
               / calc_h(r_points_, index + 1)
               * params_->v1(0);
  } else if (index == r_points_.size() - 1) {
    return params_->f(r_points_[index], 0) * calc_cross_h(r_points_, index)
            + params_->v2(0);
  } else {
    return params_->f(r_points_[index], 0)  * calc_cross_h(r_points_, index);
  }
  assert(false);
}

/* clang-format on */
