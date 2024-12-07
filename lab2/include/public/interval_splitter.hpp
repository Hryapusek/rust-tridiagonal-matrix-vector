#include <vector>
#include <cstdio>

#include <contract/contract.hpp>

#include <defines.hpp>

auto split_interval(const Number_t& left, const Number_t& right, size_t num_intervals) -> std::vector<Number_t>;

// Calculate the length of an interval `index-1` to `index`
auto calc_h(const std::vector<Number_t>& intervals, size_t index) -> Number_t;

// Calculate the cross h of an interval
auto calc_cross_h(const std::vector<Number_t>& intervals, size_t index) -> Number_t;

/// @return middle point between `index` and `index - 1`
auto middle_point(const std::vector<Number_t>& intervals, size_t index) -> Number_t;
