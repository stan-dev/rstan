#ifndef STAN_MATH_PRIM_MAT_META_IS_CONSTANT_STRUCT_HPP
#define STAN_MATH_PRIM_MAT_META_IS_CONSTANT_STRUCT_HPP

#include <stan/math/prim/arr/meta/is_constant_struct.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_constant.hpp>

namespace stan {

template <typename T, int R, int C>
struct is_constant_struct<Eigen::Matrix<T, R, C> > {
  enum { value = is_constant_struct<T>::value };
};

template <typename T>
struct is_constant_struct<Eigen::Block<T> > {
  enum { value = is_constant_struct<T>::value };
};

}  // namespace stan
#endif
