# IgGeneUsage

<details>

* Version: 1.6.0
* GitHub: https://github.com/snaketron/IgGeneUsage
* Source code: https://github.com/cran/IgGeneUsage
* Date/Publication: 2021-05-19
* Number of recursive dependencies: 99

Run `revdep_details(, "IgGeneUsage")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in ‘IgGeneUsage-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: DGU
    > ### Title: Differential gene usage in immune repertoires
    > ### Aliases: DGU
    > 
    > ### ** Examples
    > 
    > # input data
    ...
    Chain 1: Gradient evaluation took 0.000119 seconds
    Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.19 seconds.
    Chain 1: Adjust your expectations accordingly!
    Chain 1: 
    Chain 1: 
    Chain 1: Iteration:    1 / 1500 [  0%]  (Warmup)
    Chain 1: Iteration:  250 / 1500 [ 16%]  (Warmup)
    Chain 1: Iteration:  500 / 1500 [ 33%]  (Warmup)
    Chain 1: Iteration:  501 / 1500 [ 33%]  (Sampling)
    Killed
    ```

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
       2.   ├─rstan::sampling(...)
       3.   └─rstan::sampling(...)
       4.     └─rstan:::.local(object, ...)
       5.       └─base::sapply(nfits, FUN = function(x) x@mode == 0)
       6.         └─base::lapply(X = X, FUN = FUN, ...)
       7.           └─rstan:::FUN(X[[i]], ...)
      
      [ FAIL 1 | WARN 1 | SKIP 0 | PASS 92 ]
      Error: Test failures
      Execution halted
      
       *** caught segfault ***
      address 0x18, cause 'memory not mapped'
      An irrecoverable exception occurred. R is aborting now ...
      Segmentation fault (core dumped)
    ```

## In both

*   checking for hidden files and directories ... NOTE
    ```
    Found the following hidden files and directories:
      .travis.yml
    These were most likely included in error. See section ‘Package
    structure’ in the ‘Writing R Extensions’ manual.
    ```

# PosteriorBootstrap

<details>

* Version: 0.1.1
* GitHub: https://github.com/alan-turing-institute/PosteriorBootstrap
* Source code: https://github.com/cran/PosteriorBootstrap
* Date/Publication: 2021-05-14 13:02:12 UTC
* Number of recursive dependencies: 99

Run `revdep_details(, "PosteriorBootstrap")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
       1. ├─base::suppressWarnings(...) test_anpl.R:177:2
       2. │ └─base::withCallingHandlers(...)
       3. ├─rstan::vb(...)
       4. └─rstan::vb(...)
       5.   └─rstan:::.local(object, ...)
       6.     └─sampler$call_sampler(c(args, dotlist))
      
      [ FAIL 1 | WARN 0 | SKIP 3 | PASS 51 ]
      Error: Test failures
      Execution halted
      
       *** caught segfault ***
      address 0x18, cause 'memory not mapped'
      An irrecoverable exception occurred. R is aborting now ...
      Segmentation fault (core dumped)
    ```

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘Rcpp’ ‘StanHeaders’ ‘dplyr’ ‘ggplot2’ ‘gridExtra’ ‘rstan’ ‘tibble’
      All declared Imports should be used.
    ```

# ProbReco

<details>

* Version: 0.1.0.1
* GitHub: https://github.com/anastasiospanagiotelis/ProbReco
* Source code: https://github.com/cran/ProbReco
* Date/Publication: 2020-09-24 08:10:06 UTC
* Number of recursive dependencies: 67

Run `revdep_details(, "ProbReco")` for more info

</details>

## Newly broken

*   checking whether package ‘ProbReco’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/ProbReco/new/ProbReco.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking installed package size ... NOTE
    ```
      installed size is 19.1Mb
      sub-directories of 1Mb or more:
        libs  18.7Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘Rdpack’
      All declared Imports should be used.
    ```

## Installation

### Devel

```
* installing *source* package ‘ProbReco’ ...
** package ‘ProbReco’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/new/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/BH/include'    -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c RcppExports.cpp -o RcppExports.o
In file included from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/Eigen/Core:397,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/Eigen/Dense:1,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/RcppEigenForward.h:30,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/RcppEigen.h:25,
                 from RcppExports.cpp:4:
...
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
Error: package or namespace load failed for ‘ProbReco’ in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/ProbReco/new/ProbReco.Rcheck/00LOCK-ProbReco/00new/ProbReco/libs/ProbReco.so':
  /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/ProbReco/new/ProbReco.Rcheck/00LOCK-ProbReco/00new/ProbReco/libs/ProbReco.so: undefined symbol: _ZN3tbb8internal26task_scheduler_observer_v37observeEb
Error: loading failed
Execution halted
ERROR: loading failed
* removing ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/ProbReco/new/ProbReco.Rcheck/ProbReco’


```
### CRAN

```
* installing *source* package ‘ProbReco’ ...
** package ‘ProbReco’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/old/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/BH/include'    -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c RcppExports.cpp -o RcppExports.o
In file included from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/Eigen/Core:397,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/Eigen/Dense:1,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/RcppEigenForward.h:30,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include/RcppEigen.h:25,
                 from RcppExports.cpp:4:
...
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (ProbReco)


```
# prophet

<details>

* Version: 1.0
* GitHub: https://github.com/facebook/prophet
* Source code: https://github.com/cran/prophet
* Date/Publication: 2021-03-30 12:10:07 UTC
* Number of recursive dependencies: 95

Run `revdep_details(, "prophet")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      Chain 4: 
      Chain 4:  Elapsed Time: 6.521 seconds (Warm-up)
      Chain 4:                1.079 seconds (Sampling)
      Chain 4:                7.6 seconds (Total)
      Chain 4: 
      [ FAIL 0 | WARN 1 | SKIP 0 | PASS 375 ]
      > 
      > proc.time()
         user  system elapsed 
      137.126   3.012 140.123 
      
       *** caught segfault ***
      address 0x18, cause 'memory not mapped'
      An irrecoverable exception occurred. R is aborting now ...
      Segmentation fault (core dumped)
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 53.9Mb
      sub-directories of 1Mb or more:
        libs  53.1Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘StanHeaders’ ‘methods’ ‘rstantools’
      All declared Imports should be used.
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# rstanarm

<details>

* Version: 2.21.1
* GitHub: https://github.com/stan-dev/rstanarm
* Source code: https://github.com/cran/rstanarm
* Date/Publication: 2020-07-20 18:20:03 UTC
* Number of recursive dependencies: 136

Run `revdep_details(, "rstanarm")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      17: test_code(NULL, exprs, env)
      18: source_file(path, child_env(env), wrap = wrap)
      19: FUN(X[[i]], ...)
      20: lapply(test_paths, test_one_file, env = env, wrap = wrap)
      21: doTryCatch(return(expr), name, parentenv, handler)
      22: tryCatchOne(expr, names, parentenv, handlers[[1L]])
      23: tryCatchList(expr, classes, parentenv, handlers)
      24: tryCatch(code, testthat_abort_reporter = function(cnd) {    cat(conditionMessage(cnd), "\n")    NULL})
      25: with_reporter(reporters$multi, lapply(test_paths, test_one_file,     env = env, wrap = wrap))
      26: test_files(test_dir = test_dir, test_package = test_package,     test_paths = test_paths, load_helpers = load_helpers, reporter = reporter,     env = env, stop_on_failure = stop_on_failure, stop_on_warning = stop_on_warning,     wrap = wrap, load_package = load_package)
      27: test_files(test_dir = path, test_paths = test_paths, test_package = package,     reporter = reporter, load_helpers = load_helpers, env = env,     stop_on_failure = stop_on_failure, stop_on_warning = stop_on_warning,     wrap = wrap, load_package = load_package, parallel = parallel)
      28: test_dir("testthat", package = package, reporter = reporter,     ..., load_package = "installed")
      29: test_check("rstanarm", invert = FALSE, filter = if (Sys.getenv("NOT_CRAN") !=     "true") "stan_functions")
      An irrecoverable exception occurred. R is aborting now ...
      Segmentation fault (core dumped)
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 396.3Mb
      sub-directories of 1Mb or more:
        R       1.6Mb
        doc     4.5Mb
        libs  389.6Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘RcppParallel’
      All declared Imports should be used.
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# stanette

<details>

* Version: 2.21.2
* GitHub: https://github.com/stan-dev/rstan
* Source code: https://github.com/cran/stanette
* Date/Publication: 2021-08-16 08:30:02 UTC
* Number of recursive dependencies: 115

Run `revdep_details(, "stanette")` for more info

</details>

## Newly broken

*   checking whether package ‘stanette’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/stanette/new/stanette.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking installed package size ... NOTE
    ```
      installed size is 357.2Mb
      sub-directories of 1Mb or more:
        lib    65.7Mb
        libs  289.9Mb
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘stanette’ ...
** package ‘stanette’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
echo "make libdparse.a ..."
make libdparse.a ...
(cd d; make CC="gcc -std=gnu99" CFLAGS="-g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fPIC" libdparse.a)
make[1]: Entering directory '/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/stanette/new/stanette.Rcheck/00_pkg_src/stanette/src/d'
gcc -std=gnu99 -g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fPIC   -c -o arg.o arg.c
gcc -std=gnu99 -g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fPIC   -c -o parse.o parse.c
...
/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/stanette/RcppEigen/include/Eigen/src/Core/ProductEvaluators.h:124:75:   required from ‘Eigen::internal::product_evaluator<Eigen::Product<Lhs, Rhs, Option>, ProductTag, LhsShape, RhsShape>::product_evaluator(const XprType&) [with Lhs = Eigen::Product<Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<double, double>, const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, const Eigen::Matrix<double, 1, -1> >, const Eigen::Transpose<Eigen::Matrix<double, -1, 1> > >, Eigen::Matrix<double, -1, -1>, 0>; Rhs = Eigen::Matrix<double, -1, 1>; int Options = 0; int ProductTag = 6; LhsShape = Eigen::DenseShape; RhsShape = Eigen::DenseShape; typename Eigen::internal::traits<typename Eigen::Product<Lhs, Rhs, Option>::Rhs>::Scalar = double; typename Eigen::internal::traits<typename Eigen::Product<Lhs, Rhs, Option>::Lhs>::Scalar = double; Eigen::internal::product_evaluator<Eigen::Product<Lhs, Rhs, Option>, ProductTag, LhsShape, RhsShape>::XprType = Eigen::Product<Eigen::Product<Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<double, double>, const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, const Eigen::Matrix<double, 1, -1> >, const Eigen::Transpose<Eigen::Matrix<double, -1, 1> > >, Eigen::Matrix<double, -1, -1>, 0>, Eigen::Matrix<double, -1, 1>, 0>]’
/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/stanette/RcppEigen/include/Eigen/src/Core/ProductEvaluators.h:35:90:   required from ‘Eigen::internal::evaluator<Eigen::Product<Lhs, Rhs, Option> >::evaluator(const XprType&) [with Lhs = Eigen::Product<Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<double, double>, const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, const Eigen::Matrix<double, 1, -1> >, const Eigen::Transpose<Eigen::Matrix<double, -1, 1> > >, Eigen::Matrix<double, -1, -1>, 0>; Rhs = Eigen::Matrix<double, -1, 1>; int Options = 0; Eigen::internal::evaluator<Eigen::Product<Lhs, Rhs, Option> >::XprType = Eigen::Product<Eigen::Product<Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<double, double>, const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, const Eigen::Matrix<double, 1, -1> >, const Eigen::Transpose<Eigen::Matrix<double, -1, 1> > >, Eigen::Matrix<double, -1, -1>, 0>, Eigen::Matrix<double, -1, 1>, 0>]’
/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/stanette/RcppEigen/include/Eigen/src/Core/Product.h:132:22:   required from ‘Eigen::internal::dense_product_base<Lhs, Rhs, Option, 6>::operator const Scalar() const [with Lhs = Eigen::Product<Eigen::CwiseBinaryOp<Eigen::internal::scalar_product_op<double, double>, const Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, const Eigen::Matrix<double, 1, -1> >, const Eigen::Transpose<Eigen::Matrix<double, -1, 1> > >, Eigen::Matrix<double, -1, -1>, 0>; Rhs = Eigen::Matrix<double, -1, 1>; int Option = 0; Eigen::internal::dense_product_base<Lhs, Rhs, Option, 6>::Scalar = double]’
./stan/mcmc/hmc/hamiltonians/dense_e_metric.hpp:23:56:   required from ‘double stan::mcmc::dense_e_metric<Model, BaseRNG>::T(stan::mcmc::dense_e_point&) [with Model = stan::model::model_base; BaseRNG = boost::random::additive_combine_engine<boost::random::linear_congruential_engine<unsigned int, 40014, 0, 2147483563>, boost::random::linear_congruential_engine<unsigned int, 40692, 0, 2147483399> >]’
./stan/mcmc/hmc/hamiltonians/dense_e_metric.hpp:22:10:   required from here
/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/stanette/RcppEigen/include/Eigen/src/Core/DenseCoeffsBase.h:55:30: warning: ignoring attributes on template argument ‘Eigen::internal::packet_traits<double>::type’ {aka ‘__m128d’} [-Wignored-attributes]
/usr/lib/R/etc/Makeconf:177: recipe for target 'stan_fit.o' failed
make: *** [stan_fit.o] Error 1
ERROR: compilation failed for package ‘stanette’
* removing ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/stanette/new/stanette.Rcheck/stanette’


```
### CRAN

```
* installing *source* package ‘stanette’ ...
** package ‘stanette’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
echo "make libdparse.a ..."
make libdparse.a ...
(cd d; make CC="gcc -std=gnu99" CFLAGS="-g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fPIC" libdparse.a)
make[1]: Entering directory '/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/stanette/old/stanette.Rcheck/00_pkg_src/stanette/src/d'
gcc -std=gnu99 -g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fPIC   -c -o arg.o arg.c
gcc -std=gnu99 -g -O2 -fdebug-prefix-map=/build/r-base-5XUBcI/r-base-4.1.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fPIC   -c -o parse.o parse.c
...
Creating a new generic function for ‘as.mcmc.list’ in package ‘stanette’
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (stanette)


```
