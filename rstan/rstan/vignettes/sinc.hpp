template <typename T0__>
T0__
sinc(const T0__& x, std::ostream* pstream__) {
  return x != 0.0 ? sin(x) / x : (x + 1.0);
}

stan::math::vari
sinc(const stan::math::var& x, std::ostream* pstream__) {
  double x_ = x.val();
  double f = x_ != 0.0 ? sin(x_) / x_ : 1.0;
  double dfdx_ = x_ != 0.0 ? (cos(x_) - sin(x_)) / x_ : 0.0;
  return precomp_v_vari(f, x.vi_, dfdx_);
}

