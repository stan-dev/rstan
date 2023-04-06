# beanz

<details>

* Version: 2.4
* GitHub: NA
* Source code: https://github.com/cran/beanz
* Date/Publication: 2018-11-05 21:30:03 UTC
* Number of recursive dependencies: 94

Run `revdepcheck::revdep_details(, "beanz")` for more info

</details>

## Newly broken

*   checking whether package ‘beanz’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/beanz/new/beanz.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking installed package size ... NOTE
    ```
      installed size is  7.1Mb
      sub-directories of 1Mb or more:
        libs   5.6Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘rstantools’
      All declared Imports should be used.
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘beanz’ ...
** package ‘beanz’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/bs.stan
Wrote C++ file "stan_files/bs.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c stan_files/bs.cc -o stan_files/bs.o
In file included from stan_files/bs.cc:3:
...
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [stan_files/bs.o] Error 1
rm stan_files/bs.cc
ERROR: compilation failed for package ‘beanz’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/beanz/new/beanz.Rcheck/beanz’


```
### CRAN

```
* installing *source* package ‘beanz’ ...
** package ‘beanz’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/bs.stan
Wrote C++ file "stan_files/bs.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c stan_files/bs.cc -o stan_files/bs.o
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/ds.stan
...
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (beanz)


```
# dfpk

<details>

* Version: 3.5.1
* GitHub: https://github.com/artemis-toumazi/dfpk
* Source code: https://github.com/cran/dfpk
* Date/Publication: 2018-11-09 15:20:06 UTC
* Number of recursive dependencies: 50

Run `revdepcheck::revdep_details(, "dfpk")` for more info

</details>

## Newly broken

*   checking whether package ‘dfpk’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/dfpk/new/dfpk.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘dfpk’ ...
** package ‘dfpk’ successfully unpacked and MD5 sums checked
** using staged installation
** libs


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c Modules.cpp -o Modules.o
In file included from Modules.cpp:3:
In file included from ./include/models.hpp:5:
In file included from /Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/rstan/include/rstan/rstaninc.hpp:4:
...
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:66:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [Modules.o] Error 1
ERROR: compilation failed for package ‘dfpk’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/dfpk/new/dfpk.Rcheck/dfpk’


```
### CRAN

```
* installing *source* package ‘dfpk’ ...
** package ‘dfpk’ successfully unpacked and MD5 sums checked
** using staged installation
** libs


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c Modules.cpp -o Modules.o


clang -arch arm64 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -c init.c -o init.o
...
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (dfpk)


```
# idem

<details>

* Version: 5.1
* GitHub: https://github.com/olssol/idem
* Source code: https://github.com/cran/idem
* Date/Publication: 2021-01-27 09:40:02 UTC
* Number of recursive dependencies: 104

Run `revdepcheck::revdep_details(, "idem")` for more info

</details>

## Newly broken

*   checking whether package ‘idem’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/idem/new/idem.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘idem’ ...
** package ‘idem’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/idem.stan
Wrote C++ file "stan_files/idem.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c stan_files/idem.cc -o stan_files/idem.o
In file included from stan_files/idem.cc:3:
...
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [stan_files/idem.o] Error 1
rm stan_files/idem.cc
ERROR: compilation failed for package ‘idem’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/idem/new/idem.Rcheck/idem’


```
### CRAN

```
* installing *source* package ‘idem’ ...
** package ‘idem’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/idem.stan
Wrote C++ file "stan_files/idem.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c stan_files/idem.cc -o stan_files/idem.o

...
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (idem)


```
# oncomsm

<details>

* Version: 0.1.3
* GitHub: https://github.com/Boehringer-Ingelheim/oncomsm
* Source code: https://github.com/cran/oncomsm
* Date/Publication: 2023-03-11 10:20:02 UTC
* Number of recursive dependencies: 125

Run `revdepcheck::revdep_details(, "oncomsm")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      
      `actual`:   FALSE
      `expected`: TRUE 
      Backtrace:
          ▆
       1. └─oncomsm (local) test_calibration(scale_factor, shape) at test-sampling.R:221:6
       2.   └─testthat::expect_true(...) at test-sampling.R:213:4
      
      [ FAIL 2 | WARN 0 | SKIP 3 | PASS 57 ]
      Deleting unused snapshots:
      • plots/plot-mstate-srp-model-2.svg
      • plots/plot-mstate-srp-model-3.svg
      • plots/plot-srp-model-2.svg
      Error: Test failures
      Execution halted
    ```

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘rstantools’
      All declared Imports should be used.
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# visit

<details>

* Version: 2.1
* GitHub: NA
* Source code: https://github.com/cran/visit
* Date/Publication: 2019-08-26 09:30:02 UTC
* Number of recursive dependencies: 91

Run `revdepcheck::revdep_details(, "visit")` for more info

</details>

## Newly broken

*   checking whether package ‘visit’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/visit/new/visit.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘sqldf’
      All declared Imports should be used.
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘visit’ ...
** package ‘visit’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/visit.stan
Wrote C++ file "stan_files/visit.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c stan_files/visit.cc -o stan_files/visit.o
In file included from stan_files/visit.cc:3:
...
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [stan_files/visit.o] Error 1
rm stan_files/visit.cc
ERROR: compilation failed for package ‘visit’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/visit/new/visit.Rcheck/visit’


```
### CRAN

```
* installing *source* package ‘visit’ ...
** package ‘visit’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/visit.stan
Wrote C++ file "stan_files/visit.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -Wno-nonnull -c stan_files/visit.cc -o stan_files/visit.o

...
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (visit)


```
