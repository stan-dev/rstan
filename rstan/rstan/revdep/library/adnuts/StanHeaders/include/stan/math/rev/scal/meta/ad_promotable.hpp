#ifndef STAN_MATH_REV_SCAL_META_AD_PROMOTABLE_HPP
#define STAN_MATH_REV_SCAL_META_AD_PROMOTABLE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/meta/ad_promotable.hpp>

namespace stan {
namespace math {

template <typename T>
struct ad_promotable<T, var> {
  enum { value = ad_promotable<T, double>::value };
};

template <>
struct ad_promotable<var, var> {
  enum { value = true };
};

}  // namespace math
}  // namespace stan
#endif
