#pragma once

#include <functional>

// First argument is r, second is t
using Function_type = std::function<double(double, double)>;

using T_function_type = std::function<double(double)>;
using R_function_type = std::function<double(double)>;

using Number_t = double;

struct InputParameters {
  Number_t Rl;
  Number_t Rr; // [Rl, Rr]
  Number_t T; // [0, T]

  // First type condition
  T_function_type v1; // u(rL) = v1(t)

  // Third type condition
  Number_t hi2; // -k * du/dr = hi2 * u(rR) - v2(t)
  R_function_type phi;
  T_function_type v2;

  // Just input functions
  Function_type k;
  Function_type q;
  Function_type f;
};
