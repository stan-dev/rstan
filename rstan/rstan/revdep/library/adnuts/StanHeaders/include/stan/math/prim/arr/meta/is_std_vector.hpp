#ifndef STAN_MATH_PRIM_ARR_META_IS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_ARR_META_IS_STD_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_std_vector.hpp>
#include <vector>

namespace stan {
template <typename T>
struct is_std_vector<std::vector<T> > {
  enum { value = 1 };
  typedef T type;
};
}  // namespace stan

#endif
