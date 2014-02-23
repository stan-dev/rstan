#include <gtest/gtest.h>
#include <rstan/stan_fit.hpp>
#include <stan/common/do_print.hpp>


TEST(do_print, valid) {
  for (int n = -10; n <= 10; n++) 
    for (int refresh = -10; refresh <= 10; refresh++)
      for (int last = -10; last <= 10; last++) {
        EXPECT_EQ(rstan::do_print(n, refresh, last),
                  stan::common::do_print(n, n==last, refresh));
      }
}
