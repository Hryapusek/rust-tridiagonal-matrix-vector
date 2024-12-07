#pragma once

#include <memory>

#include <interface/i_main_matrix_calculator.hpp>
#include <input_parameters.hpp>

class DefaultMainMatrixCalculator : public IMainMatrixCalculator
{
 public:
  auto calc_a(size_t index) -> Number_t override;
  auto calc_b(size_t index) -> Number_t override;
  auto calc_c(size_t index) -> Number_t override;
  auto calc_g(size_t index) -> Number_t override;

 protected:
  std::vector<Number_t> intervals_;
  std::shared_ptr<InputParameters> params_;
};
