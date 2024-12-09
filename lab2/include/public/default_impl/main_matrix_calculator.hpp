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
    std::vector<Number_t> r_points
  )
    : params_(params)
    , r_points_(std::move(r_points))
  {}

  auto calc_a(size_t index) const -> Number_t override;
  auto calc_b(size_t index) const -> Number_t override;
  auto calc_c(size_t index) const -> Number_t override;
  auto calc_g(size_t index) const -> Number_t override;

  auto r_points() const -> std::vector<Number_t> const& { return r_points_; }

  auto params() const -> std::shared_ptr<InputParameters> const& { return params_; }

 protected:
  std::vector<Number_t> r_points_;
  std::shared_ptr<InputParameters> params_;
};
