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
    precondition(index < intervals_.size(), "index out of range");
  };
  return middle_point(intervals_, index) * params_->k(middle_point(intervals_, index), 0)
       / calc_cross_h(intervals_, index) / intervals_.at(index) / calc_h(intervals_, index);
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
      index != intervals_.size() - 1,
      "You should not calculate b for index == intervals_.size() - 1"
    );
    precondition(index < intervals_.size(), "index out of range");
  };
  auto up = middle_point(intervals_, index + 1)
          * params_->k(middle_point(intervals_, index + 1), 0);
  auto down = calc_cross_h(intervals_, index) * intervals_.at(index)
            * calc_h(intervals_, index + 1);
  return up / down;
}

/* clang-format off */
auto DefaultMainMatrixCalculator::calc_c(size_t index) const -> Number_t
{
  contract(fun) {
    precondition(index != 0, "You should not calculate anything for index == 0 - you already have v function");
    precondition(index < intervals_.size(), "index out of range");
  };
  if (index == intervals_.size() - 1) {
    return - middle_point(intervals_, index) * params_->k(middle_point(intervals_, index), 0)
                 / calc_h(intervals_, index)
                 / calc_cross_h(intervals_, index)
                 / intervals_[index]
            - params_->hi2
                 / calc_cross_h(intervals_, index)
            - params_->q(intervals_[index], 0);
  } else {
    return -middle_point(intervals_, index) * params_->k(middle_point(intervals_, index), 0)
                 / calc_h(intervals_, index + 1)
                 / calc_cross_h(intervals_, index)
                 / intervals_[index]
            - middle_point(intervals_, index + 1) * params_->k(middle_point(intervals_, index + 1), 0)
                 / calc_h(intervals_, index + 1)
                 / calc_cross_h(intervals_, index)
                 / intervals_[index]
            - params_->q(intervals_[index], 0);
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
    return params_->f(intervals_[index], 0) +
               middle_point(intervals_, index) * params_->k(middle_point(intervals_, index), 0)
               / calc_cross_h(intervals_, index) 
               / intervals_[index] 
               / calc_h(intervals_, index + 1);
  } else if (index == intervals_.size() - 1) {
    return params_->f(intervals_[index], 0)
            + params_->v2(0) / calc_cross_h(intervals_, index);
  } else {
    return params_->f(intervals_[index], 0);
  }
  assert(false);
}

/* clang-format on */
