double
sinc(const double& x, std::ostream* pstream__) {
  return x != 0.0 ? sin(x) / x : 1.0;
}

var
sinc(const var& x, std::ostream* pstream__) {
  double x_ = x.val();
  double f = x_ != 0.0 ? sin(x_) / x_ : 1.0;
  double dfdx_ = x_ != 0.0 ? (cos(x_) - sin(x_)) / x_ : 0.0;
  return precomputed_gradients(f, {x.vi_}, {dfdx_});
}

