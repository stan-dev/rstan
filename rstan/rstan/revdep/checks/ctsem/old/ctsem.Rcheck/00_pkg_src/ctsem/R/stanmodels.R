# Part of the ctsem package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# This file is only intended to be used during the installation process
# nocov start
MODELS_HOME <- "src"
if (!file.exists(MODELS_HOME)) MODELS_HOME <- sub("R$", "src", getwd())

w32 <-  .Machine$sizeof.pointer == 4 #.Platform$OS.type == "windows" &&
# w32 <- FALSE
# 
# # stan_files <- file.path(MODELS_HOME, 
# #   ifelse(w32, "stan_files/ctsmW32.stan", "stan_files/ctsm.stan" ))
# # 
# # if(!file.exists(paste0('./src/stan_files/ctsm',ifelse(w32,'W32',''),'.stan')) {
# #   file.rename(paste0('./src/stan_files/ctsm',ifelse(w32,'W32',''),'.bak'), 
# #   paste0('./src/stan_files/ctsm',ifelse(w32,'W32',''),'.stan')) #rename backup .stan file to use
# # }
# 
# # suppressMessages(try(file.rename('./src/stan_files/ctsm.stan','./src/stan_files/ctsm.del'),silent = TRUE))
# file.copy(paste0('./src/stan_files/ctsm',ifelse(w32,'W32',''),'.prep'),
#   paste0('./src/stan_files/ctsm',ifelse(w32,'W32',''),'.stan'),overwrite=TRUE,copy.date=TRUE) #rename unused .stan file to avoid errors

stan_files <- dir(file.path(MODELS_HOME, paste0('stan_files',ifelse(w32,'32',''))),
                  pattern = "stan$", full.names = TRUE)
stanmodels <- lapply(stan_files, function(f) {
  model_cppname <- sub("\\.stan$", "", basename(f))
  stanfit <- suppressWarnings(rstan::stanc(f, allow_undefined = TRUE, 
                          obfuscate_model_name = FALSE))
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
                            model_cppcode = stanfit$cppcode)
  return(do.call(methods::new, args = c(stanfit[-(1:3)], Class = "stanmodel", 
                 mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))
  }
)
names(stanmodels) <- sub("\\.stan$", "", basename(stan_files))
rm(MODELS_HOME)
rm(w32)
# nocov end
