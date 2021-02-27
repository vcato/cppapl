#include <iostream>
#include "vector.hpp"


template <typename T>
inline std::ostream& operator<<(std::ostream& s, const vector<T> &v)
{
  s << "[";

  for (const auto &x : v) {
    s << " " << x;
  }

  s << " ]";
  return s;
}
