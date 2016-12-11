# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015 Trustees of Columbia University
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

rstan_ess <- function(sim, n) {
  # Args:
  #   n: Chain index starting from 1.
  ess <- .Call(effective_sample_size, sim, n - 1)
  ess
} 

rstan_splitrhat <- function(sim, n) {
  # Args:
  #   n: Chain index starting from 1.
  if (sim$n_save[1] - sim$warmup2[1] < 2) return(NaN)
  rhat <- .Call(split_potential_scale_reduction, sim, n - 1)
  rhat
}

rstan_splitrhat2_cpp <- function(sims) {
  # Args:
  #   sim: samples of several chains _without_ warmup
  rhat <- .Call(split_potential_scale_reduction2, sims)
  rhat
}

rstan_ess2_cpp <- function(sims) {
  # Args:
  #   sim: samples of several chains _without_ warmup
  ess <- .Call(effective_sample_size2, sims)
  ess
} 

rstan_seq_perm <- function(n, chains, seed, chain_id = 1) {
  # Args:
  #   n: length of sequence to be generated 
  #   chains: the number of chains, for which the permuations are applied
  #   seed: the seed for RNG 
  #   chain_id: the chain id, for which the returned permuation is applied 
  # 
  conf <- list(n = n, chains = chains, seed = seed, chain_id = chain_id) 
  perm <- .Call(seq_permutation, conf)
  perm + 1L # start from 1 
} 
