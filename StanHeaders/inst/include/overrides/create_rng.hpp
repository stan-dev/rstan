#ifndef STAN_SERVICES_UTIL_CREATE_RNG_HPP
#define STAN_SERVICES_UTIL_CREATE_RNG_HPP

#include <boost/random/mixmax.hpp>
#include <boost/random/additive_combine.hpp>

namespace stan {

#ifndef NEW_RSTAN
using rng_t = boost::ecuyer1988;
#else
using rng_t = boost::random::mixmax;
#endif

namespace services {
namespace util {

#ifndef NEW_RSTAN
/**
 * Creates a pseudo random number generator from a random seed
 * and a chain id by initializing the PRNG with the seed and
 * then advancing past pow(2, 50) times the chain ID draws to
 * ensure different chains sample from different segments of the
 * pseudo random number sequence.
 *
 * Chain IDs should be kept to larger values than one to ensure
 * that the draws used to initialized transformed data are not
 * duplicated.
 *
 * @param[in] seed the random seed
 * @param[in] chain the chain id
 * @return a boost::ecuyer1988 instance
 */
inline boost::ecuyer1988 create_rng(unsigned int seed, unsigned int chain) {
  using boost::uintmax_t;
  static constexpr uintmax_t DISCARD_STRIDE = static_cast<uintmax_t>(1) << 50;
  boost::ecuyer1988 rng(seed);
  // always discard at least 1 to avoid issue with small seeds for certain RNG
  // distributions. See stan#3167 and boostorg/random#92
  rng.discard(std::max(static_cast<uintmax_t>(1), DISCARD_STRIDE * chain));
  return rng;
}
#else
/**
 * Creates a pseudo random number generator from a random seed
 * and a chain id by initializing the PRNG with the seed and
 * then advancing past pow(2, 50) times the chain ID draws to
 * ensure different chains sample from different segments of the
 * pseudo random number sequence.
 *
 * Chain IDs should be kept to larger values than one to ensure
 * that the draws used to initialized transformed data are not
 * duplicated.
 *
 * @param[in] seed the random seed
 * @param[in] chain the chain id
 * @return an stan::rng_t instance
 */
inline rng_t create_rng(unsigned int seed, unsigned int chain) {
  // RNG state is 128 bits, but user only provides 64 total bits
  // Additionally, there are issues if all 128 bits are 0, hence
  // the 1 as the second argument
  rng_t rng(0, 1, seed, chain);
  return rng;
}
#endif

}  // namespace util
}  // namespace services
}  // namespace stan
#endif
