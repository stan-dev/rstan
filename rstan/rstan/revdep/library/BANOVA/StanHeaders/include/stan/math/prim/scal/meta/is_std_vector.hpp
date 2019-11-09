#ifndef STAN_MATH_PRIM_SCAL_META_IS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_STD_VECTOR_HPP

#include <vector>

namespace stan {
template <typename T>
struct is_std_vector {
  enum { value = 0 };
  typedef T type;
};
}  // namespace stan

#endif
