#ifndef STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_DIAG_HPP
#define STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_DIAG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_vector.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, Eigen::Dynamic>
quad_form_diag(const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& mat,
               const Eigen::Matrix<T2, R, C>& vec) {
  using boost::math::tools::promote_args;
  check_vector("quad_form_diag", "vec", vec);
  check_square("quad_form_diag", "mat", mat);
  int size = vec.size();
  check_size_match("quad_form_diag", "rows of mat", mat.rows(), "size of vec",
                   size);
  Eigen::Matrix<typename promote_args<T1, T2>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      result(size, size);
  for (int i = 0; i < size; i++) {
    result(i, i) = vec(i) * vec(i) * mat(i, i);
    for (int j = i + 1; j < size; ++j) {
      typename promote_args<T1, T2>::type temp = vec(i) * vec(j);
      result(j, i) = temp * mat(j, i);
      result(i, j) = temp * mat(i, j);
    }
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
