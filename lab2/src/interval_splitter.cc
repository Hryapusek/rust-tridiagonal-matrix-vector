#include <interval_splitter.hpp>

auto split_interval(Number_t const& left, Number_t const& right, size_t num_intervals) -> std::vector<Number_t>
{
  std::vector<Number_t> intervals;

  contract(fun)
  {
    precondition(num_intervals > 0, "invalid number of intervals");
  };

  auto interval_size = (right - left) / num_intervals;
  for(size_t i = 0; i < num_intervals; ++i) {
    intervals.push_back(left + interval_size * i);
  }
  intervals.push_back(right);
  return intervals;
}

// Calculate the length of an interval `index-1` to `index`
auto calc_h(std::vector<Number_t> const& points, size_t index) -> Number_t
{
  contract(fun)
  {
    precondition(index < points.size(), "index out of range");
    precondition(index > 0, "h can not be calculated for the first interval");
  };

  return (points.at(index) - points.at(index - 1));
}

// Calculate the cross h of an interval
auto calc_cross_h(std::vector<Number_t> const& points, size_t index) -> Number_t
{
  if(index == 0) {
    return calc_h(points, 1) / 2;
  }
  else if(index == points.size() - 1) {
    return calc_h(points, index) / 2;
  }
  else {
    return (calc_h(points, index) + calc_h(points, index + 1)) / 2;
  }
}

/// @return middle point between `index` and `index - 1`
auto middle_point(std::vector<Number_t> const& points, size_t index) -> Number_t
{
  return (points.at(index) + points.at(index - 1)) / 2;
}
