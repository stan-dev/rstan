#ifndef STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/make_nu.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/lkj_corr_lpdf.hpp>

namespace stan {
namespace math {

// LKJ_Corr(L|eta) [ L Cholesky factor of correlation matrix
//                  eta > 0; eta == 1 <-> uniform]
template <bool propto, typename T_covar, typename T_shape>
inline return_type_t<T_covar, T_shape> lkj_corr_cholesky_lpdf(
    const T_covar& L, const T_shape& eta) {
  const unsigned int K = L.rows();
  if (K == 0) {
    return 0.0;
  }
  using lp_ret = return_type_t<T_covar, T_shape>;
  static constexpr const char* function = "lkj_corr_cholesky_lpdf";
  check_positive(function, "Shape parameter", eta);

  const auto& L_ref = to_ref(L);
  check_cholesky_factor(function, "Random variable", L_ref);


  lp_ret lp(0.0);

  if constexpr (include_summand<propto, T_shape>::value) {
    lp += do_lkj_constant(eta, K);
  }
  if constexpr (include_summand<propto, T_covar, T_shape>::value) {
    const int Km1 = K - 1;
    Eigen::Matrix<value_type_t<T_covar>, Eigen::Dynamic, 1> log_diagonals
        = log(L_ref.diagonal().tail(Km1).array());
    Eigen::Matrix<lp_ret, Eigen::Dynamic, 1> values(Km1);
    for (int k = 0; k < Km1; k++) {
      values(k) = (Km1 - k - 1) * log_diagonals(k);
    }
    if constexpr (is_constant_all<scalar_type<T_shape>>::value) {
      if (eta == 1.0) {
        lp += sum(values);
        return (lp);
      }
    }
    values += multiply(2.0 * eta - 2.0, log_diagonals);
    lp += sum(values);
  }

  return lp;
}

template <typename T_covar, typename T_shape>
inline return_type_t<T_covar, T_shape> lkj_corr_cholesky_lpdf(
    const T_covar& L, const T_shape& eta) {
  return lkj_corr_cholesky_lpdf<false>(L, eta);
}

}  // namespace math
}  // namespace stan
#endif
