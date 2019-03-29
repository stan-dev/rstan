#ifndef BOOST_RANDOM_R_HPP
#define BOOST_RANDOM_R_HPP
#include <R.h>
#include <cstdint>

class boost_random_R {
public:
  boost_random_R() = default;
  typedef uint32_t result_type;
  static const bool has_fixed_range = true;  
  result_type min() {return 1;};
  result_type max() {return 2147483562;};
  result_type operator()() {
    GetRNGstate();
    result_type out = static_cast<result_type>(min() + 
      unif_rand() * (max() - min()));
    PutRNGstate();
    return out;
  };
  void discard(unsigned long long j) {
    GetRNGstate();
    for (unsigned long long k = 0; k < j; k++)
      double discarded = unif_rand();
    PutRNGstate();
  };
  
  template<typename T_>
  void seed(T_) {}; // not allowed to touch R's PRNG state except from R
  
  template<typename T0_, typename T1_>
  friend bool operator==(T0_, T1_) {return false;};
  template<typename T0_, typename T1_>
  friend bool operator!=(T0_, T1_) {return true;};
};
#endif
