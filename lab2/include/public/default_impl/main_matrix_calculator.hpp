#pragma once

#include <memory>

#include <Eigen/Dense>

#include <interface/i_main_matrix_calculator.hpp>
#include <input_parameters.hpp>

class DefaultMainMatrixCalculator : public IMainMatrixCalculator
{
 public:
  explicit DefaultMainMatrixCalculator(
    std::shared_ptr<InputParameters> params,
    std::vector<Number_t> intervals
  )
    : params_(params)
    , intervals_(std::move(intervals))
  {}

  auto calc_a(size_t index) const -> Number_t override;
  auto calc_b(size_t index) const -> Number_t override;
  auto calc_c(size_t index) const -> Number_t override;
  auto calc_g(size_t index) const -> Number_t override;

  auto intervals() const -> std::vector<Number_t> const& { return intervals_; }

  auto params() const -> std::shared_ptr<InputParameters> const& { return params_; }

 protected:
  std::vector<Number_t> intervals_;
  std::shared_ptr<InputParameters> params_;
};
