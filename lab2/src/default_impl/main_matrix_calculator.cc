#include <default_impl/main_matrix_calculator.hpp>

#include <cassert>

#include <contract/contract.hpp>

#include <interval_splitter.hpp>

auto DefaultMainMatrixCalculator::calc_a(size_t r_index, Number_t t) const -> Number_t
{
  contract(fun)
  {
    precondition(
      r_index != 0,
      "You should not calculate anything for index == 0 - you already have v function"
    );
    precondition(r_index != 1, "You should not calculate a for index == 1");
    precondition(r_index < r_points_.size(), "index out of range");
  };
  auto mid_point = middle_point(r_points_, r_index);
  auto k_value = params_->k(mid_point, t);
  auto h_value = calc_h(r_points_, r_index);

  return mid_point * k_value / h_value;
}

auto DefaultMainMatrixCalculator::calc_b(size_t r_index, Number_t t) const -> Number_t
{
  contract(fun)
  {
    precondition(
      r_index != 0,
      "You should not calculate anything for index == 0 - you already have v function"
    );
    precondition(
      r_index != r_points_.size() - 1,
      "You should not calculate b for index == intervals_.size() - 1"
    );
    precondition(r_index < r_points_.size(), "index out of range");
  };
  auto up = middle_point(r_points_, r_index + 1) * params_->k(middle_point(r_points_, r_index + 1), t);
  auto down = calc_h(r_points_, r_index + 1);
  return up / down;
}

/* clang-format off */
auto DefaultMainMatrixCalculator::calc_c(size_t r_index, Number_t t) const -> Number_t
{
  contract(fun) {
    precondition(r_index != 0, "You should not calculate anything for index == 0 - you already have v function");
    precondition(r_index < r_points_.size(), "index out of range");
  };
  if (r_index == r_points_.size() - 1) {
    return - middle_point(r_points_, r_index) * params_->k(middle_point(r_points_, r_index), t)
                 / calc_h(r_points_, r_index)
            - params_->hi2
            - params_->q(r_points_[r_index], t)
                 * calc_cross_h(r_points_, r_index);
  } else {
    auto m_point_index = middle_point(r_points_, r_index);
    auto m_point_index_plus_1 = middle_point(r_points_, r_index + 1);
    auto k_index = params_->k(m_point_index, t);
    auto k_index_plus_1 = params_->k(m_point_index_plus_1, t);
    auto h_index_plus_1 = calc_h(r_points_, r_index + 1);
    auto cross_h_index = calc_cross_h(r_points_, r_index);
    auto r = r_points_[r_index];
    auto q_index = params_->q(r_index, t);

    return -m_point_index * k_index / h_index_plus_1
           -m_point_index_plus_1 * k_index_plus_1 / h_index_plus_1
           -q_index * calc_cross_h(r_points_, r);
  }
  assert(false);
}

/* clang-format on */

/* clang-format off */
auto DefaultMainMatrixCalculator::calc_g(size_t r_index, Number_t t) const -> Number_t {
  contract(fun) {
    precondition(r_index != 0, "You should not calculate anything for index == 0 - you already have v function");
  };
  if (r_index == 1) {
    return params_->f(r_points_[r_index], t) +
               middle_point(r_points_, r_index) 
               * params_->k(middle_point(r_points_, r_index), t)
               / calc_h(r_points_, r_index + 1)
               * params_->v1(0);
  } else if (r_index == r_points_.size() - 1) {
    return params_->f(r_points_[r_index], t) * calc_cross_h(r_points_, r_index)
            + params_->v2(0);
  } else {
    return params_->f(r_points_[r_index], t)  * calc_cross_h(r_points_, r_index);
  }
  assert(false);
}

/* clang-format on */
