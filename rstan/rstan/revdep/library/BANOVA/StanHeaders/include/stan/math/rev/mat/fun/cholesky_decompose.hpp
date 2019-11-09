#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <algorithm>

namespace stan {
namespace math {

class cholesky_block : public vari {
 public:
  int M_;
  int block_size_;
  typedef Eigen::Block<Eigen::MatrixXd> Block_;
  vari** variRefA_;
  vari** variRefL_;

  /**
   * Constructor for cholesky function.
   *
   * Stores varis for A.  Instantiates and stores varis for L.
   * Instantiates and stores dummy vari for upper triangular part of var
   * result returned in cholesky_decompose function call
   *
   * variRefL aren't on the chainable autodiff stack, only used for storage
   * and computation. Note that varis for L are constructed externally in
   * cholesky_decompose.
   *
   * block_size_ determined using the same calculation Eigen/LLT.h
   *
   * @param A matrix
   * @param L_A matrix, cholesky factor of A
   */
  cholesky_block(const Eigen::Matrix<var, -1, -1>& A,
                 const Eigen::Matrix<double, -1, -1>& L_A)
      : vari(0.0),
        M_(A.rows()),
        variRefA_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            A.rows() * (A.rows() + 1) / 2)),
        variRefL_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            A.rows() * (A.rows() + 1) / 2)) {
    size_t pos = 0;
    block_size_ = std::max((M_ / 8 / 16) * 16, 8);
    block_size_ = std::min(block_size_, 128);
    for (size_type j = 0; j < M_; ++j) {
      for (size_type i = j; i < M_; ++i) {
        variRefA_[pos] = A.coeffRef(i, j).vi_;
        variRefL_[pos] = new vari(L_A.coeffRef(i, j), false);
        ++pos;
      }
    }
  }

  /**
   * Symbolic adjoint calculation for cholesky factor A
   *
   * @param L cholesky factor
   * @param Lbar matrix of adjoints of L
   */
  inline void symbolic_rev(Block_& L, Block_& Lbar) {
    using Eigen::Lower;
    using Eigen::StrictlyUpper;
    using Eigen::Upper;
    L.transposeInPlace();
    Lbar = (L * Lbar.triangularView<Lower>()).eval();
    Lbar.triangularView<StrictlyUpper>()
        = Lbar.adjoint().triangularView<StrictlyUpper>();
    L.triangularView<Upper>().solveInPlace(Lbar);
    L.triangularView<Upper>().solveInPlace(Lbar.transpose());
  }

  /**
   * Reverse mode differentiation algorithm refernce:
   *
   * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
   *
   */
  virtual void chain() {
    using Eigen::Block;
    using Eigen::Lower;
    using Eigen::MatrixXd;
    using Eigen::StrictlyUpper;
    using Eigen::Upper;
    MatrixXd Lbar(M_, M_);
    MatrixXd L(M_, M_);

    Lbar.setZero();
    L.setZero();
    size_t pos = 0;
    for (size_type j = 0; j < M_; ++j) {
      for (size_type i = j; i < M_; ++i) {
        Lbar.coeffRef(i, j) = variRefL_[pos]->adj_;
        L.coeffRef(i, j) = variRefL_[pos]->val_;
        ++pos;
      }
    }

    for (int k = M_; k > 0; k -= block_size_) {
      int j = std::max(0, k - block_size_);
      Block_ R = L.block(j, 0, k - j, j);
      Block_ D = L.block(j, j, k - j, k - j);
      Block_ B = L.block(k, 0, M_ - k, j);
      Block_ C = L.block(k, j, M_ - k, k - j);
      Block_ Rbar = Lbar.block(j, 0, k - j, j);
      Block_ Dbar = Lbar.block(j, j, k - j, k - j);
      Block_ Bbar = Lbar.block(k, 0, M_ - k, j);
      Block_ Cbar = Lbar.block(k, j, M_ - k, k - j);
      if (Cbar.size() > 0) {
        Cbar = D.transpose()
                   .triangularView<Upper>()
                   .solve(Cbar.transpose())
                   .transpose();
        Bbar.noalias() -= Cbar * R;
        Dbar.noalias() -= Cbar.transpose() * C;
      }
      symbolic_rev(D, Dbar);
      Rbar.noalias() -= Cbar.transpose() * B;
      Rbar.noalias() -= Dbar.selfadjointView<Lower>() * R;
      Dbar.diagonal() *= 0.5;
      Dbar.triangularView<StrictlyUpper>().setZero();
    }
    pos = 0;
    for (size_type j = 0; j < M_; ++j)
      for (size_type i = j; i < M_; ++i)
        variRefA_[pos++]->adj_ += Lbar.coeffRef(i, j);
  }
};

class cholesky_scalar : public vari {
 public:
  int M_;
  vari** variRefA_;
  vari** variRefL_;

  /**
   * Constructor for cholesky function.
   *
   * Stores varis for A Instantiates and stores varis for L Instantiates
   * and stores dummy vari for upper triangular part of var result returned
   * in cholesky_decompose function call
   *
   * variRefL aren't on the chainable autodiff stack, only used for storage
   * and computation. Note that varis for L are constructed externally in
   * cholesky_decompose.
   *
   * @param A matrix
   * @param L_A matrix, cholesky factor of A
   */
  cholesky_scalar(const Eigen::Matrix<var, -1, -1>& A,
                  const Eigen::Matrix<double, -1, -1>& L_A)
      : vari(0.0),
        M_(A.rows()),
        variRefA_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            A.rows() * (A.rows() + 1) / 2)),
        variRefL_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            A.rows() * (A.rows() + 1) / 2)) {
    size_t accum = 0;
    size_t accum_i = accum;
    for (size_type j = 0; j < M_; ++j) {
      for (size_type i = j; i < M_; ++i) {
        accum_i += i;
        size_t pos = j + accum_i;
        variRefA_[pos] = A.coeffRef(i, j).vi_;
        variRefL_[pos] = new vari(L_A.coeffRef(i, j), false);
      }
      accum += j;
      accum_i = accum;
    }
  }

  /**
   * Reverse mode differentiation algorithm refernce:
   *
   * Mike Giles. An extended collection of matrix derivative results for
   * forward and reverse mode AD.  Jan. 2008.
   *
   * Note algorithm  as laid out in Giles is row-major, so Eigen::Matrices
   * are explicitly storage order RowMajor, whereas Eigen defaults to
   * ColumnMajor. Also note algorithm starts by calculating the adjoint for
   * A(M_ - 1, M_ - 1), hence pos on line 94 is decremented to start at pos
   * = M_ * (M_ + 1) / 2.
   */
  virtual void chain() {
    using Eigen::Matrix;
    using Eigen::RowMajor;
    Matrix<double, -1, -1, RowMajor> adjL(M_, M_);
    Matrix<double, -1, -1, RowMajor> LA(M_, M_);
    Matrix<double, -1, -1, RowMajor> adjA(M_, M_);
    size_t pos = 0;
    for (size_type i = 0; i < M_; ++i) {
      for (size_type j = 0; j <= i; ++j) {
        adjL.coeffRef(i, j) = variRefL_[pos]->adj_;
        LA.coeffRef(i, j) = variRefL_[pos]->val_;
        ++pos;
      }
    }

    --pos;
    for (int i = M_ - 1; i >= 0; --i) {
      for (int j = i; j >= 0; --j) {
        if (i == j) {
          adjA.coeffRef(i, j) = 0.5 * adjL.coeff(i, j) / LA.coeff(i, j);
        } else {
          adjA.coeffRef(i, j) = adjL.coeff(i, j) / LA.coeff(j, j);
          adjL.coeffRef(j, j)
              -= adjL.coeff(i, j) * LA.coeff(i, j) / LA.coeff(j, j);
        }
        for (int k = j - 1; k >= 0; --k) {
          adjL.coeffRef(i, k) -= adjA.coeff(i, j) * LA.coeff(j, k);
          adjL.coeffRef(j, k) -= adjA.coeff(i, j) * LA.coeff(i, k);
        }
        variRefA_[pos--]->adj_ += adjA.coeffRef(i, j);
      }
    }
  }
};

/**
 * Reverse mode specialization of cholesky decomposition
 *
 * Internally calls Eigen::LLT rather than using
 * stan::math::cholesky_decompose in order to use an inplace decomposition.
 *
 * Note chainable stack varis are created below in Matrix<var, -1, -1>
 *
 * @param A Matrix
 * @return L cholesky factor of A
 */
inline Eigen::Matrix<var, -1, -1> cholesky_decompose(
    const Eigen::Matrix<var, -1, -1>& A) {
  check_square("cholesky_decompose", "A", A);
  check_symmetric("cholesky_decompose", "A", A);

  Eigen::Matrix<double, -1, -1> L_A(value_of_rec(A));
  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower> L_factor(L_A);
  check_pos_definite("cholesky_decompose", "m", L_factor);

  // Memory allocated in arena.
  // cholesky_scalar gradient faster for small matrices compared to
  // cholesky_block
  vari* dummy = new vari(0.0, false);
  Eigen::Matrix<var, -1, -1> L(A.rows(), A.cols());
  if (L_A.rows() <= 35) {
    cholesky_scalar* baseVari = new cholesky_scalar(A, L_A);
    size_t accum = 0;
    size_t accum_i = accum;
    for (size_type j = 0; j < L.cols(); ++j) {
      for (size_type i = j; i < L.cols(); ++i) {
        accum_i += i;
        size_t pos = j + accum_i;
        L.coeffRef(i, j).vi_ = baseVari->variRefL_[pos];
      }
      for (size_type k = 0; k < j; ++k)
        L.coeffRef(k, j).vi_ = dummy;
      accum += j;
      accum_i = accum;
    }
  } else {
    cholesky_block* baseVari = new cholesky_block(A, L_A);
    size_t pos = 0;
    for (size_type j = 0; j < L.cols(); ++j) {
      for (size_type i = j; i < L.cols(); ++i) {
        L.coeffRef(i, j).vi_ = baseVari->variRefL_[pos++];
      }
      for (size_type k = 0; k < j; ++k)
        L.coeffRef(k, j).vi_ = dummy;
    }
  }
  return L;
}
}  // namespace math
}  // namespace stan
#endif
