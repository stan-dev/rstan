template <typename T0__, typename T1__>
typename boost::math::tools::promote_args<T0__, T1__>::type
besselK(const T0__& v, const T1__& z, std::ostream* pstream__) {
  return boost::math::cyl_bessel_k(v, z);
}
