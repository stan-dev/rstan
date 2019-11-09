#ifndef STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/meta/scalar_type.hpp>

namespace stan {

template <typename T, int R, int C>
struct scalar_type<Eigen::Matrix<T, R, C> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T, int R, int C>
struct scalar_type<const Eigen::Matrix<T, R, C> > {
  typedef typename scalar_type<T>::type type;
};

template <typename T, int R, int C>
struct scalar_type<Eigen::Matrix<T, R, C>&> {
  typedef typename scalar_type<T>::type type;
};

template <typename T, int R, int C>
struct scalar_type<const Eigen::Matrix<T, R, C>&> {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
struct scalar_type<Eigen::Block<T> > {
  typedef typename scalar_type<T>::type type;
};
}  // namespace stan
#endif
