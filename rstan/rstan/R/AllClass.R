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

setClass(Class = "cxxdso",
         representation = representation(
           sig = "list", # A list of function signature that would be returned by cxxfuncion 
           dso_saved = "logical", # flag for if the dso is saved or not
           # dso_last_path = 'character', # where the dso is saved last time 
           dso_filename = "character", # the dso file name when it is created the first time
           modulename = 'character', # in rstan, we always compile the c++ code with an Rcpp module
           system = "character", # what is the OS (R.version$system)?  
           cxxflags = "character", # the CXXFLAGS used to compile the DSO 
           .CXXDSOMISC = "environment" # an envir to save 
                                 #  1. the returned object by cxxfuncion using name cxxfun 
                                 #  2. the file path used last time using name dso_last_path 
                                 #  3. The binary dso with name dso_bin, which is a raw vector.  
                                 #     We put it here since the environment is not copied 
                                 #     when assigned to another 
                                 #     http://cran.r-project.org/doc/manuals/R-lang.html#Environment-objects
         ),
         validity = function(object) {
           length(object@sig) > 0 && identical(object@system, R.version$system)
         })

setClass(Class = "stanmodel",
         representation = representation(
           model_name = "character",
           model_code = "character",
           model_cpp = "list", 
             # model_cppname (used to define Rcpp module)  & 
             # model_cppcode (just the C++ code for the model) 
           mk_cppmodule = "function", # use a function to create the cpp module
           dso = 'cxxdso'), 
         validity = function(object) {
           return(TRUE)
         })

setClass(Class = "stanmodel_with_data", 
         representation = representation(model = "ANY", # FIXME
                                         data = "list"),
         validity = function(object) return(TRUE) )

setMethod("initialize", "stanmodel_with_data", 
          function(.Object, model, data) {
            
            stan_fit_cpp_module <- model@mk_cppmodule(model)
            cxxfun <- rstan:::grab_cxxfun(model@dso)
            smwd <- try(new(stan_fit_cpp_module, data))
            if (is(smwd, "try-error")) {
              stop("something went horribly wrong")
            }
            .Object@model <- smwd
            .Object@data <- data
            return(.Object)
})
         
setClass(Class = "stanfit",
         representation = representation(
           model_name = "character", 
           model_pars = "character", 
           par_dims = "list", # the order of par in the list should match those in model_pars
           mode = "integer", # 0: samples; 1: test_grad (no samples); 2: other error (no samples) 
           sim = "list", # samples with other helper information: pars_oi, dims_oi, fnames_oi, n_flatnames,
           inits = "list", 
           stan_args = "list", 
           stanmodel = "stanmodel", # the instance of S4 class stanmodel 
           date = "character", # the date samples were drawn 
           .MISC = "environment"
         ),  
         validity = function(object) {
           if(length(object@sim) > 0 && !is.null(object@sim$samples)) {
             NAs <- rapply(object@sim$samples, f = function(x) any(is.na(x)))
             if(any(NAs)) return(paste("The following variables have undefined values: ",
                                       unique(names(NAs[NAs])), collapse = ","))
           }
           return(TRUE) 
         })

# list all methods for stanfit 
# showMethods(class = 'stanfit', where = getNamespace("rstan"))
