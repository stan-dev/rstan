double sinc(const double& x, std::ostream* pstream__) {
  return x != 0.0 ? sin(x) / x : 1.0;
}

stan::math::var sinc(const stan::math::var& x, std::ostream* pstream__) {
  return stan::math::make_callback_vari(sinc(x.val(), pstream__), [x](const auto& vi) mutable {
    double dfdx_ = x.val() != 0.0 ? (cos(x.val()) - sin(x.val())) / x.val() : 0.0;
    x.adj() += vi.adj() * dfdx_;
  });
}

