#ifndef STAN_MATH_PRIM_ARR_META_CONTAINS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_ARR_META_CONTAINS_STD_VECTOR_HPP

#include <stan/math/prim/arr/meta/is_std_vector.hpp>

namespace stan {

template <typename T1, typename T2 = double, typename T3 = double,
          typename T4 = double, typename T5 = double, typename T6 = double>
struct contains_std_vector {
  enum {
    value = is_std_vector<T1>::value || is_std_vector<T2>::value
            || is_std_vector<T3>::value || is_std_vector<T4>::value
            || is_std_vector<T5>::value || is_std_vector<T6>::value
  };
};

}  // namespace stan
#endif
