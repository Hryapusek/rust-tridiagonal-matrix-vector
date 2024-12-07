#pragma once

#include <cstdio>

using Number_t = double;

class IMainMatrixCalculator
{
 public:
  using Number_t = double;
  virtual auto calc_a(size_t index) -> Number_t = 0;
  virtual auto calc_b(size_t index) -> Number_t = 0;
  virtual auto calc_c(size_t index) -> Number_t = 0;
  virtual auto calc_g(size_t index) -> Number_t = 0;
};
