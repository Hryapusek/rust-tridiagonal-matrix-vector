#pragma once

#include <defines.hpp>

struct InputParameters {
  Number_t Rl;
  Number_t Rr; // [Rl, Rr]
  Number_t T; // [0, T]

  // First type condition
  T_Function_type v1; // u(rL) = v1(t)

  // Third type condition
  Number_t hi2; // -k * du/dr = hi2 * u(rR) - v2(t)
  R_Function_type phi;
  T_Function_type v2;

  // Just input functions
  R_T_Function_type k;
  R_T_Function_type q;
  R_T_Function_type f;
};
