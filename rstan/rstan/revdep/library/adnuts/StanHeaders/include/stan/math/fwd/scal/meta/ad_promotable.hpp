#ifndef STAN_MATH_FWD_SCAL_META_AD_PROMOTABLE_HPP
#define STAN_MATH_FWD_SCAL_META_AD_PROMOTABLE_HPP

#include <stan/math/prim/scal/meta/ad_promotable.hpp>

namespace stan {
namespace math {

template <typename T>
struct fvar;

template <typename V, typename T>
struct ad_promotable<V, fvar<T> > {
  enum { value = ad_promotable<V, T>::value };
};

}  // namespace math
}  // namespace stan
#endif
