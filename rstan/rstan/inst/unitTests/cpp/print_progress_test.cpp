#include <gtest/gtest.h>
#include <rstan/stan_fit.hpp>
#include <stan/common/print_progress.hpp>
#include <sstream>


TEST(print_progress, old) {
  std::stringstream rstan;
  std::stringstream stan;
  rstan::print_progress(0, 100, 1, true, rstan);
  stan::common::print_progress(0, 0, 100, 1, true, "\r", "", stan);
  EXPECT_EQ(stan.str(), rstan.str());

  for (int m = 0; m < 100; m++) {
    rstan.str("");
    stan.str("");
    rstan::print_progress(0, 100, 1, true, rstan);
    stan::common::print_progress(0, 0, 100, 1, true, "\r", "", stan);
    EXPECT_EQ(stan.str(), rstan.str());
  }

  for (int m = 0; m < 100; m++) {
    rstan.str("");
    stan.str("");
    rstan::print_progress(0, 100, 1, false, rstan);
    stan::common::print_progress(0, 0, 100, 1, false, "\r", "", stan);
    EXPECT_EQ(stan.str(), rstan.str());
  }

  for (int m = 0; m < 100; m++) {
    rstan.str("");
    stan.str("");
    rstan::print_progress(0, 100, 5, false, rstan);
    stan::common::print_progress(0, 0, 100, 5, false, "\r", "", stan);
    EXPECT_EQ(stan.str(), rstan.str());
  }

}
