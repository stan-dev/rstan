# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2023 Trustees of Columbia University
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

library(dplyr)
library(stringr)

# created via: ./stanc --dump-stan-math-signatures
signatures <- file.path("tools", "signatures.txt")
functions <- file.path("tools", "stan-functions.txt")

system2("sed", args = paste("'s/((.*=> vector,/(vector_function,/g'", signatures),
        stdout = functions, stderr = NULL)
stopifnot(.Last.value == 0)
system2("sed", args = paste("-i 's/((.*=> real,/(real_function,/g'", functions),
        stdout = NULL, stderr = NULL)
stopifnot(.Last.value == 0)
system2("sed", args = paste("-i 's/((.*=> array\\[\\] real,/(real_array_function,/g'", functions),
        stdout = NULL, stderr = NULL)
stopifnot(.Last.value == 0)
system2("sed", args = paste("-i 's/=> tuple(/=> tuple{/g'", functions),
        stdout = NULL, stderr = NULL)
stopifnot(.Last.value == 0)
system2("sed", args = paste("-i 's/ => /;/g'", functions),
        stdout = NULL, stderr = NULL)
stopifnot(.Last.value == 0)
system2("sed", args = paste("-i 's/(/;(/g'", functions),
        stdout = NULL, stderr = NULL)
stopifnot(.Last.value == 0)
system2("sed", args = paste("-i '1iStanFunction;Arguments;ReturnType'", functions),
        stdout = NULL, stderr = NULL)
stopifnot(.Last.value == 0)
system2("sed", args = paste("-i '1i# This file is semicolon delimited'", functions),
        stdout = NULL, stderr = NULL)
stopifnot(.Last.value == 0)

rosetta <- cbind(RFunction = NA_character_,
                 read.table(functions, header = TRUE, sep = ";", 
                            quote = NULL, stringsAsFactors = FALSE, strip.white = TRUE)) %>% 
  mutate(RFunction = if_else(StanFunction %in% unlist(sapply(search(), ls)),
                             StanFunction, NA_character_),
         templated = str_detect(Arguments, fixed(ReturnType)),
         Arguments = case_when(templated & paste0("(", ReturnType, ")") == Arguments ~ "(T)",
                               templated ~ str_replace(Arguments, 
                                                       fixed(paste0(", ", ReturnType, ")")),
                                                       ", T)"),
                               TRUE ~ Arguments),
         ReturnType = if_else(str_ends(Arguments, fixed("T)")), "T", ReturnType)) %>% 
  group_by(StanFunction, templated) %>% 
  slice_head(n = 1) %>% 
  ungroup %>% 
  arrange(str_to_lower(StanFunction), desc(templated)) %>% 
  group_by(StanFunction, RFunction) %>% 
  summarize(Arguments = paste(Arguments, collapse = ";"), 
            ReturnType = paste(ReturnType, collapse = ";"),
            .groups = "drop") %>% 
  arrange(str_to_lower(StanFunction)) %>% 
  as.data.frame

rosetta <- within(rosetta, RFunction[StanFunction == "add"] <- "+")  
rosetta <- within(rosetta, RFunction[grepl("algebra_solve", StanFunction)] <- "uniroot")
rosetta <- within(rosetta, RFunction[StanFunction == "append_array"] <- "c")
rosetta <- within(rosetta, RFunction[StanFunction == "append_col"] <- "cbind")  
rosetta <- within(rosetta, RFunction[StanFunction == "append_row"] <- "rbind")
rosetta <- within(rosetta, RFunction[StanFunction == "arg"] <- "Arg")  
rosetta <- within(rosetta, RFunction[grepl("^bernoulli_.*cdf", StanFunction)] <- "pbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "bernoulli"] <- "dbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "bernoulli_logit"] <- "dbinom")
rosetta <- within(rosetta, RFunction[grepl("^bernoulli_.*pmf", StanFunction)] <- "dbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "bernoulli_rng"] <- "rbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "bessel_first_kind"] <- "besselJ")
rosetta <- within(rosetta, RFunction[StanFunction == "bessel_second_kind"] <- "besselY")
rosetta <- within(rosetta, RFunction[StanFunction == "beta"] <- "dbeta")
rosetta <- within(rosetta, RFunction[grepl("^beta_[lc]*cdf", StanFunction)] <- "pbeta")
rosetta <- within(rosetta, RFunction[StanFunction == "beta_lpdf"] <- "dbeta")
rosetta <- within(rosetta, RFunction[StanFunction == "beta_rng"] <- "rbeta")
rosetta <- within(rosetta, RFunction[grepl("^binomial_[lc]*cdf", StanFunction)] <- "pbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "binomial_coefficient_log"] <- "lchoose")
rosetta <- within(rosetta, RFunction[grepl("^binomial_logit", StanFunction)] <- "binomial")
rosetta <- within(rosetta, RFunction[StanFunction == "binomial"] <- "dbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "binomial_lpmf"] <- "dbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "binomial_rng"] <- "rbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "block"] <- "subset")
rosetta <- within(rosetta, RFunction[grepl("^categorical", StanFunction)] <- "dmultinom")
rosetta <- within(rosetta, RFunction[StanFunction == "categorical_rng"] <- "rmultinom")
rosetta <- within(rosetta, RFunction[grepl("^cauchy[_lpdf]*", StanFunction)] <- "dcauchy")
rosetta <- within(rosetta, RFunction[grepl("^cauchy_[lc]*cdf", StanFunction)] <- "pcauchy")
rosetta <- within(rosetta, RFunction[StanFunction == "cauchy_rng"] <- "rcauchy")
rosetta <- within(rosetta, RFunction[StanFunction == "ceil"] <- "ceiling")
rosetta <- within(rosetta, RFunction[grepl("^chi_square[_lpdf]*", StanFunction)] <- "dchisq")
rosetta <- within(rosetta, RFunction[grepl("^chi_square_[lc]*cdf", StanFunction)] <- "pchisq")
rosetta <- within(rosetta, RFunction[StanFunction == "chi_square_rng"] <- "rchisq")
rosetta <- within(rosetta, RFunction[StanFunction == "cholesky_decompose"] <- "chol")
rosetta <- within(rosetta, RFunction[StanFunction == "choose"] <- "choose")
rosetta <- within(rosetta, RFunction[StanFunction == "col"] <- "subset")
rosetta <- within(rosetta, RFunction[StanFunction == "cols"] <- "NCOL")
rosetta <- within(rosetta, RFunction[grepl("^columns_", StanFunction)] <- "apply")
rosetta <- within(rosetta, RFunction[StanFunction == "conj"] <- "Conj")
rosetta <- within(rosetta, RFunction[StanFunction == "cumulative_sum"] <- "cumsum")
rosetta <- within(rosetta, RFunction[StanFunction == "determinant"] <- "det")
rosetta <- within(rosetta, RFunction[StanFunction == "diag_matrix"] <- "diag")
rosetta <- within(rosetta, RFunction[StanFunction == "diagonal"] <- "diag")
rosetta <- within(rosetta, RFunction[StanFunction == "dims"] <- "dim")
rosetta <- within(rosetta, RFunction[StanFunction == "distance"] <- "dist")
rosetta <- within(rosetta, RFunction[StanFunction == "divide"] <- "/")
rosetta <- within(rosetta, RFunction[StanFunction == "dot_self"] <- "crossprod")
rosetta <- within(rosetta, RFunction[StanFunction == "dot_product"] <- "%*%")
rosetta <- within(rosetta, RFunction[StanFunction == "e"] <- "exp")
rosetta <- within(rosetta, RFunction[grepl("^eigen", StanFunction)] <- "eigen")
rosetta <- within(rosetta, RFunction[grepl("^erf", StanFunction)] <- "pnorm")
rosetta <- within(rosetta, RFunction[grepl("^exponential[_lpdf]*", StanFunction)] <- "dexp")
rosetta <- within(rosetta, RFunction[grepl("^exponential_[lc]*cdf", StanFunction)] <- "pexp")
rosetta <- within(rosetta, RFunction[StanFunction == "exponential_rng"] <- "rexp")
rosetta <- within(rosetta, RFunction[StanFunction == "fabs"] <- "abs")
rosetta <- within(rosetta, RFunction[StanFunction == "fmax"] <- "max")
rosetta <- within(rosetta, RFunction[StanFunction == "fmin"] <- "min")
rosetta <- within(rosetta, RFunction[StanFunction == "fmod"] <- "%%")
rosetta <- within(rosetta, RFunction[grepl("^gamma[_lpdf]*", StanFunction)] <- "dgamma")
rosetta <- within(rosetta, RFunction[grepl("^gamma_[lc]*cdf", StanFunction)] <- "pgamma")
rosetta <- within(rosetta, RFunction[StanFunction == "gamma_rng"] <- "rgamma")
rosetta <- within(rosetta, RFunction[grepl("^gamma_[pq]$", StanFunction)] <- "pgamma")
rosetta <- within(rosetta, RFunction[StanFunction == "generalized_inverse"] <- "MASS::ginv")
rosetta <- within(rosetta, RFunction[grepl("^hypergeometric[_lpmf]*", StanFunction)] <- "dhyper")
rosetta <- within(rosetta, RFunction[grepl("^hypergeometric_[lc]*cdf", StanFunction)] <- "phyper")
rosetta <- within(rosetta, RFunction[StanFunction == "hypergeometric_rng"] <- "rhyper")
rosetta <- within(rosetta, RFunction[StanFunction == "if_else"] <- "ifelse")
rosetta <- within(rosetta, RFunction[StanFunction == "integrate_1d"] <- "integrate")
rosetta <- within(rosetta, RFunction[StanFunction == "inverse"] <- "solve")
rosetta <- within(rosetta, RFunction[StanFunction == "inverse_spd"] <- "solve")
rosetta <- within(rosetta, RFunction[StanFunction == "inv_logit"] <- "plogis")
rosetta <- within(rosetta, RFunction[StanFunction == "inv_Phi"] <- "qnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "is_inf"] <- "is.finite")
rosetta <- within(rosetta, RFunction[StanFunction == "is_nan"] <- "is.nan")
rosetta <- within(rosetta, RFunction[StanFunction == "lgamma"] <- "lgamma")
rosetta <- within(rosetta, RFunction[StanFunction == "log1m"] <- "log1p")
rosetta <- within(rosetta, RFunction[StanFunction == "log_determinant"] <- "determinant")
rosetta <- within(rosetta, RFunction[StanFunction == "log_inv_logit"] <- "plogis")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_and"] <- "&&")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_eq"] <- "==")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_gt"] <- ">")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_gte"] <- ">=")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_lt"] <- "<")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_lte"] <- "<=")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_negation"] <- "!")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_neq"] <- "!=")
rosetta <- within(rosetta, RFunction[StanFunction == "logical_or"] <- "||")
rosetta <- within(rosetta, RFunction[grepl("^logistic[_lpdf]*", StanFunction)] <- "dlogis")
rosetta <- within(rosetta, RFunction[grepl("^logistic_[lc]*cdf", StanFunction)] <- "plogis")
rosetta <- within(rosetta, RFunction[StanFunction == "logistic_rng"] <- "rlogis")
rosetta <- within(rosetta, RFunction[StanFunction == "logit"] <- "qlogis")
rosetta <- within(rosetta, RFunction[grepl("^lognormal[_lpdf]*", StanFunction)] <- "dlnorm")
rosetta <- within(rosetta, RFunction[grepl("^lognormal_[lc]*cdf", StanFunction)] <- "plnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "lognormal_rng"] <- "rlnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "machine_precision"] <- ".Machine")
rosetta <- within(rosetta, RFunction[StanFunction == "map_rect"] <- "vapply")
rosetta <- within(rosetta, RFunction[StanFunction == "mdivide_left"] <- "solve")
rosetta <- within(rosetta, RFunction[StanFunction == "mdivide_left_spd"] <- "solve")
rosetta <- within(rosetta, RFunction[grepl("^mdivide_left_tri_low", StanFunction)] <- "forwardsolve")
rosetta <- within(rosetta, RFunction[StanFunction == "minus"] <- "-")
rosetta <- within(rosetta, RFunction[StanFunction == "modified_bessel_first_kind"] <- "besselI")
rosetta <- within(rosetta, RFunction[StanFunction == "modified_bessel_second_kind"] <- "besselK")
rosetta <- within(rosetta, RFunction[StanFunction == "modulus"] <- "Mod")
rosetta <- within(rosetta, RFunction[grepl("^multinomial[_lpmf]*", StanFunction)] <- "dmultinom")
rosetta <- within(rosetta, RFunction[grepl("^multinomial_[lc]*cdf", StanFunction)] <- "pmultinom")
rosetta <- within(rosetta, RFunction[StanFunction == "multinomial_rng"] <- "rmultinom")
rosetta <- within(rosetta, RFunction[grepl("^multi_normal", StanFunction)] <- "mvtnorm::dmvnorm")
rosetta <- within(rosetta, RFunction[grepl("^multi_normal.*_rng$", StanFunction)] <- "mvtnorm::rmvnorm")
rosetta <- within(rosetta, RFunction[grepl("^multi_student", StanFunction)] <- "mvtnorm::dmvt")
rosetta <- within(rosetta, RFunction[grepl("^multi_student.*_rng$", StanFunction)] <- "mvtnorm::rmvt")
rosetta <- within(rosetta, RFunction[StanFunction == "multiply_lower_tri_self_transpose"] <- "crossprod")
rosetta <- within(rosetta, RFunction[StanFunction == "negative_infinity"] <- "Inf")
rosetta <- within(rosetta, RFunction[grepl("^neg_binomial[_lpmf]*", StanFunction)] <- "dnbinom")
rosetta <- within(rosetta, RFunction[grepl("^neg_binomial_[lc]*cdf", StanFunction)] <- "pnbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "neg_binomial_rng"] <- "rnbinom")
rosetta <- within(rosetta, RFunction[grepl("^neg_binomial_2[_lpmf]*", StanFunction)] <- "dnbinom")
rosetta <- within(rosetta, RFunction[grepl("^neg_binomial_2_[lc]*cdf", StanFunction)] <- "pnbinom")
rosetta <- within(rosetta, RFunction[StanFunction == "neg_binomial_2_rng"] <- "rnbinom")
rosetta <- within(rosetta, RFunction[grepl("^normal[_lpdf]*", StanFunction)] <- "dnorm")
rosetta <- within(rosetta, RFunction[grepl("^normal_[lc]*cdf", StanFunction)] <- "pnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "normal_rng"] <- "rnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "not_a_number"] <- "NaN")
rosetta <- within(rosetta, RFunction[StanFunction == "num_elements"] <- "length")
rosetta <- within(rosetta, RFunction[StanFunction == "Phi"] <- "pnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "Phi_approx"] <- "pnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "plus"] <- "+")
rosetta <- within(rosetta, RFunction[grepl("^poisson[_lpmf]*", StanFunction)] <- "dpois")
rosetta <- within(rosetta, RFunction[grepl("^poisson_[lc]*cdf", StanFunction)] <- "ppois")
rosetta <- within(rosetta, RFunction[StanFunction == "poisson_rng"] <- "rpois")
rosetta <- within(rosetta, RFunction[StanFunction == "positive_infinity"] <- "Inf")
rosetta <- within(rosetta, RFunction[StanFunction == "pow"] <- "^")
rosetta <- within(rosetta, RFunction[StanFunction == "qr_thin_Q"] <- "qr.Q")
rosetta <- within(rosetta, RFunction[StanFunction == "qr_thin_R"] <- "qr.R")
rosetta <- within(rosetta, RFunction[StanFunction == "rep_array"] <- "array")
rosetta <- within(rosetta, RFunction[StanFunction == "rep_matrix"] <- "matrix")
rosetta <- within(rosetta, RFunction[grepl("^rep_.*vector", StanFunction)] <- "rep")
rosetta <- within(rosetta, RFunction[StanFunction == "reverse"] <- "rev")
rosetta <- within(rosetta, RFunction[StanFunction == "row"] <- "subset")
rosetta <- within(rosetta, RFunction[grepl("^rows_", StanFunction)] <- "apply")
rosetta <- within(rosetta, RFunction[StanFunction == "rows"] <- "NROW")
rosetta <- within(rosetta, RFunction[StanFunction == "segment"] <- "subset")
rosetta <- within(rosetta, RFunction[StanFunction == "singular_values"] <- "svd")
rosetta <- within(rosetta, RFunction[StanFunction == "size"] <- "dim")
rosetta <- within(rosetta, RFunction[grepl("^sort_", StanFunction)] <- "sort")
rosetta <- within(rosetta, RFunction[StanFunction == "square"] <- "^")
rosetta <- within(rosetta, RFunction[StanFunction == "squared_distance"] <- "dist")
rosetta <- within(rosetta, RFunction[grepl("^std_normal_*cdf", StanFunction)] <- "pnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "std_normal_lpdf"] <- "dnorm")
rosetta <- within(rosetta, RFunction[grepl("^std_normal_*qf", StanFunction)] <- "qnorm")
rosetta <- within(rosetta, RFunction[StanFunction == "step"] <- NA_character_)
rosetta <- within(rosetta, RFunction[grepl("^student_t[_lpdf]*", StanFunction)] <- "dt")
rosetta <- within(rosetta, RFunction[grepl("^student_t_[lc]*cdf", StanFunction)] <- "pt")
rosetta <- within(rosetta, RFunction[StanFunction == "student_t_rng"] <- "rt")
rosetta <- within(rosetta, RFunction[grepl("^std_normal_*cdf", StanFunction)] <- "pnorm")
rosetta <- within(rosetta, RFunction[grepl("^svd_", StanFunction)] <- "svd")
rosetta <- within(rosetta, RFunction[StanFunction == "tgamma"] <- "gamma")
rosetta <- within(rosetta, RFunction[StanFunction == "to_array_1d"] <- "as.vector")
rosetta <- within(rosetta, RFunction[StanFunction == "to_complex"] <- "as.complex")
rosetta <- within(rosetta, RFunction[StanFunction == "to_int"] <- "as.integer")
rosetta <- within(rosetta, RFunction[StanFunction == "to_matrix"] <- "as.matrix")
rosetta <- within(rosetta, RFunction[StanFunction == "to_row_vector"] <- "as.vector")
rosetta <- within(rosetta, RFunction[StanFunction == "to_vector"] <- "as.vector")
rosetta <- within(rosetta, RFunction[StanFunction == "trace"] <- NA_character_)
rosetta <- within(rosetta, RFunction[StanFunction == "transpose"] <- "t")
rosetta <- within(rosetta, RFunction[grepl("^uniform[_lpdf]*", StanFunction)] <- "dunif")
rosetta <- within(rosetta, RFunction[grepl("^uniform_[lc]*cdf", StanFunction)] <- "punif")
rosetta <- within(rosetta, RFunction[StanFunction == "uniform_rng"] <- "runif")
rosetta <- within(rosetta, RFunction[StanFunction == "variance"] <- "var")
rosetta <- within(rosetta, RFunction[grepl("^weibull[_lpdf]*", StanFunction)] <- "dweibull")
rosetta <- within(rosetta, RFunction[grepl("^weibull_[lc]*cdf", StanFunction)] <- "pweibull")
rosetta <- within(rosetta, RFunction[StanFunction == "weibull_rng"] <- "rweibull")
rosetta <- within(rosetta, RFunction[StanFunction == "weibull_rng"] <- "rweibull")
rosetta <- within(rosetta, RFunction[StanFunction == "wishart_rng"] <- "rWishart")

save(rosetta, file = "R/sysdata.rda")
tools::resaveRdaFiles("R/sysdata.rda")
