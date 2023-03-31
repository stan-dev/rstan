double sinc(const double& x, std::ostream* pstream__) {
  return x != 0.0 ? sin(x) / x : 1.0;
}

stan::math::var sinc(const stan::math::var& x, std::ostream* pstream__) {
  double x_ = x.val();
  double f = x_ != 0.0 ? sin(x_) / x_ : 1.0;
  double dfdx_ = x_ != 0.0 ? (cos(x_) - sin(x_)) / x_ : 0.0;
  return stan::math::make_callback_vari(f, [x, dfdx_](const auto& vi) mutable {
    x.adj() += vi.adj() * dfdx_;
  });
}

