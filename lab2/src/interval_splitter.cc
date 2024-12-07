#include <interval_splitter.hpp>

auto split_interval(Number_t const& left, Number_t const& right, size_t num_intervals) -> std::vector<Number_t>
{
  std::vector<Number_t> intervals;

  contract(fun)
  {
    precondition(num_intervals > 0, "invalid number of intervals");
    postcondition(intervals.size() == num_intervals + 1, "incorrect number of intervals");
  };

  auto interval_size = (right - left) / num_intervals;
  for(size_t i = 0; i < num_intervals; ++i) {
    intervals.push_back(left + interval_size * i);
  }
  intervals.push_back(right);
  return intervals;
}

// Calculate the length of an interval `index-1` to `index`
auto calc_h(std::vector<Number_t> const& intervals, size_t index) -> Number_t
{
  contract(fun)
  {
    precondition(index < intervals.size(), "index out of range");
    precondition(index > 0, "h can not be calculated for the first interval");
  };

  return (intervals.at(index) - intervals.at(index - 1));
}

// Calculate the cross h of an interval
auto calc_cross_h(std::vector<Number_t> const& intervals, size_t index) -> Number_t
{
  if(index == 0) {
    return calc_h(intervals, 1) / 2;
  }
  else if(index == intervals.size() - 1) {
    return calc_h(intervals, index) / 2;
  }
  else {
    return (calc_h(intervals, index) + calc_h(intervals, index + 1)) / 2;
  }
}

/// @return middle point between `index` and `index - 1`
auto middle_point(std::vector<Number_t> const& intervals, size_t index) -> Number_t
{
  return (intervals.at(index) + intervals.at(index - 1)) / 2;
}
