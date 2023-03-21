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
      installed size is  7.4Mb
      sub-directories of 1Mb or more:
        libs   5.9Mb
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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/bs.cc -o stan_files/bs.o
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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/beanz/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/bs.cc -o stan_files/bs.o
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
# densEstBayes

<details>

* Version: 1.0-2.1
* GitHub: NA
* Source code: https://github.com/cran/densEstBayes
* Date/Publication: 2022-04-05 09:19:03 UTC
* Number of recursive dependencies: 50

Run `revdepcheck::revdep_details(, "densEstBayes")` for more info

</details>

## Newly broken

*   checking whether package ‘densEstBayes’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/densEstBayes/new/densEstBayes.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘methods’ ‘rstantools’
      All declared Imports should be used.
    ```

## Installation

### Devel

```
* installing *source* package ‘densEstBayes’ ...
** package ‘densEstBayes’ successfully unpacked and MD5 sums checked
** using staged installation
Registered S3 methods overwritten by 'RcppEigen':
  method               from         
  predict.fastLm       RcppArmadillo
  print.fastLm         RcppArmadillo
  summary.fastLm       RcppArmadillo
  print.summary.fastLm RcppArmadillo
Warning messages:
...
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:66:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [stanExports_PoissonSimpleMixedModel.o] Error 1
ERROR: compilation failed for package ‘densEstBayes’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/densEstBayes/new/densEstBayes.Rcheck/densEstBayes’


```
### CRAN

```
* installing *source* package ‘densEstBayes’ ...
** package ‘densEstBayes’ successfully unpacked and MD5 sums checked
** using staged installation
Registered S3 methods overwritten by 'RcppEigen':
  method               from         
  predict.fastLm       RcppArmadillo
  print.fastLm         RcppArmadillo
  summary.fastLm       RcppArmadillo
  print.summary.fastLm RcppArmadillo
Warning messages:
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
* DONE (densEstBayes)


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

*   checking installed package size ... NOTE
    ```
      installed size is  5.1Mb
      sub-directories of 1Mb or more:
        libs   4.4Mb
    ```

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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c Modules.cpp -o Modules.o
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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/dfpk/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c Modules.cpp -o Modules.o


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
# glmmPen

<details>

* Version: 1.5.3.0
* GitHub: NA
* Source code: https://github.com/cran/glmmPen
* Date/Publication: 2023-03-15 14:50:07 UTC
* Number of recursive dependencies: 93

Run `revdepcheck::revdep_details(, "glmmPen")` for more info

</details>

## Newly broken

*   checking whether package ‘glmmPen’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/glmmPen/new/glmmPen.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘glmmPen’ ...
** package ‘glmmPen’ successfully unpacked and MD5 sums checked
** using staged installation
Registered S3 methods overwritten by 'RcppEigen':
  method               from         
  predict.fastLm       RcppArmadillo
  print.fastLm         RcppArmadillo
  summary.fastLm       RcppArmadillo
  print.summary.fastLm RcppArmadillo
Warning messages:
...
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:66:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [stanExports_binomial_logit_model.o] Error 1
ERROR: compilation failed for package ‘glmmPen’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/glmmPen/new/glmmPen.Rcheck/glmmPen’


```
### CRAN

```
* installing *source* package ‘glmmPen’ ...
** package ‘glmmPen’ successfully unpacked and MD5 sums checked
** using staged installation
Registered S3 methods overwritten by 'RcppEigen':
  method               from         
  predict.fastLm       RcppArmadillo
  print.fastLm         RcppArmadillo
  summary.fastLm       RcppArmadillo
  print.summary.fastLm RcppArmadillo
Warning messages:
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
* DONE (glmmPen)


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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/idem.cc -o stan_files/idem.o
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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/idem/RcppEigen/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/idem.cc -o stan_files/idem.o

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
# nlmixr2est

<details>

* Version: 2.1.3
* GitHub: https://github.com/nlmixr2/nlmixr2est
* Source code: https://github.com/cran/nlmixr2est
* Date/Publication: 2022-11-10 17:00:19 UTC
* Number of recursive dependencies: 202

Run `revdepcheck::revdep_details(, "nlmixr2est")` for more info

</details>

## In both

*   checking whether package ‘nlmixr2est’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/nlmixr2est/new/nlmixr2est.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘nlmixr2est’ ...
** package ‘nlmixr2est’ successfully unpacked and MD5 sums checked
** using staged installation
--------[begin src/Makevars]--------
# -*- mode: makefile-gmake -*-
ARMA=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/RcppArmadillo/include
BH=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/BH/include
RCPP=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/Rcpp/include
EG=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/RcppEigen/include
SH=-isystem'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src' 
...
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:66:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [ode_cmt1.o] Error 1
ERROR: compilation failed for package ‘nlmixr2est’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/nlmixr2est/new/nlmixr2est.Rcheck/nlmixr2est’


```
### CRAN

```
* installing *source* package ‘nlmixr2est’ ...
** package ‘nlmixr2est’ successfully unpacked and MD5 sums checked
** using staged installation
--------[begin src/Makevars]--------
# -*- mode: makefile-gmake -*-
ARMA=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/RcppArmadillo/include
BH=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/BH/include
RCPP=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/Rcpp/include
EG=/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/nlmixr2est/RcppEigen/include
SH=-isystem'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src' 
...
** R
** inst
** byte-compile and prepare package for lazy loading
Error: .onLoad failed in loadNamespace() for 'rxode2', details:
  call: NULL
  error: rxode2 compiled with rxode2parse with a different solving structure
can try: install.packages('rxode2', type='source')
Execution halted
ERROR: lazy loading failed for package ‘nlmixr2est’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/nlmixr2est/old/nlmixr2est.Rcheck/nlmixr2est’


```
# OpenMx

<details>

* Version: 2.21.1
* GitHub: https://github.com/OpenMx/OpenMx
* Source code: https://github.com/cran/OpenMx
* Date/Publication: 2023-01-19 23:50:02 UTC
* Number of recursive dependencies: 146

Run `revdepcheck::revdep_details(, "OpenMx")` for more info

</details>

## Newly broken

*   checking whether package ‘OpenMx’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/OpenMx/new/OpenMx.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking installed package size ... NOTE
    ```
      installed size is 15.8Mb
      sub-directories of 1Mb or more:
        R        3.1Mb
        data     1.4Mb
        libs     4.2Mb
        models   4.7Mb
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘OpenMx’ ...
** package ‘OpenMx’ successfully unpacked and MD5 sums checked
** using staged installation
NOTE: ./configure is not an autoconf generated script.
Change default C/C++ compiler and default compile flags by editing ~/.R/Makevars
** libs
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/RcppParallel/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/rpf/include' -I/opt/R/arm64/include         -DSTRICT_R_HEADERS -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c Compute.cpp -o Compute.o
In file included from Compute.cpp:25:
In file included from ./glue.h:23:
In file included from ./omxState.h:48:
...
class omxFitFunction {
^
./omxDefines.h:131:9: note: did you mean class here?
typedef struct omxFitFunction omxFitFunction;
        ^~~~~~
        class
2 warnings and 2 errors generated.
make: *** [omxMLFitFunction.o] Error 1
ERROR: compilation failed for package ‘OpenMx’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/OpenMx/new/OpenMx.Rcheck/OpenMx’


```
### CRAN

```
* installing *source* package ‘OpenMx’ ...
** package ‘OpenMx’ successfully unpacked and MD5 sums checked
** using staged installation
NOTE: ./configure is not an autoconf generated script.
Change default C/C++ compiler and default compile flags by editing ~/.R/Makevars
** libs
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/RcppParallel/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/OpenMx/rpf/include' -I/opt/R/arm64/include         -DSTRICT_R_HEADERS -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c Compute.cpp -o Compute.o
In file included from Compute.cpp:25:
In file included from ./glue.h:23:
In file included from ./omxState.h:48:
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
* DONE (OpenMx)


```
# ProbReco

<details>

* Version: 0.1.0.1
* GitHub: https://github.com/anastasiospanagiotelis/ProbReco
* Source code: https://github.com/cran/ProbReco
* Date/Publication: 2020-09-24 08:10:06 UTC
* Number of recursive dependencies: 74

Run `revdepcheck::revdep_details(, "ProbReco")` for more info

</details>

## Newly broken

*   checking whether package ‘ProbReco’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/ProbReco/new/ProbReco.Rcheck/00install.out’ for details.
    ```

## Newly fixed

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
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/BH/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c RcppExports.cpp -o RcppExports.o
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/BH/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c scores.cpp -o scores.o
In file included from scores.cpp:7:
In file included from /Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math.hpp:19:
In file included from /Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/rev.hpp:10:
In file included from /Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/rev/fun.hpp:7:
...
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:66:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [scores.o] Error 1
ERROR: compilation failed for package ‘ProbReco’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/ProbReco/new/ProbReco.Rcheck/ProbReco’


```
### CRAN

```
* installing *source* package ‘ProbReco’ ...
** package ‘ProbReco’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/BH/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c RcppExports.cpp -o RcppExports.o
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ProbReco/BH/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c scores.cpp -o scores.o
clang++ -arch arm64 -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -L/opt/homebrew/opt/libomp/lib -o ProbReco.so RcppExports.o scores.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
ld: warning: -undefined dynamic_lookup may not work with chained fixups
installing to /Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/ProbReco/old/ProbReco.Rcheck/00LOCK-ProbReco/00new/ProbReco/libs
** R
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
# publipha

<details>

* Version: 0.1.1
* GitHub: NA
* Source code: https://github.com/cran/publipha
* Date/Publication: 2020-01-15 00:20:07 UTC
* Number of recursive dependencies: 98

Run `revdepcheck::revdep_details(, "publipha")` for more info

</details>

## Newly broken

*   checking whether package ‘publipha’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/publipha/new/publipha.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 35 marked UTF-8 strings
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘publipha’ ...
** package ‘publipha’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/cma.stan
Wrote C++ file "stan_files/cma.cc"
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"`"/Library/Frameworks/R.framework/Resources/bin/Rscript" --vanilla -e "cat(system.file('include', 'src', package = 'StanHeaders'))"`" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/cma.cc -o stan_files/cma.o
In file included from stan_files/cma.cc:3:
In file included from stan_files/cma.hpp:18:
In file included from /Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/rstan/include/rstan/rstaninc.hpp:4:
...
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [stan_files/cma.o] Error 1
rm stan_files/cma.cc
ERROR: compilation failed for package ‘publipha’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/publipha/new/publipha.Rcheck/publipha’


```
### CRAN

```
* installing *source* package ‘publipha’ ...
** package ‘publipha’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/cma.stan
Wrote C++ file "stan_files/cma.cc"
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"`"/Library/Frameworks/R.framework/Resources/bin/Rscript" --vanilla -e "cat(system.file('include', 'src', package = 'StanHeaders'))"`" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/cma.cc -o stan_files/cma.o
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/phma.stan
Wrote C++ file "stan_files/phma.cc"
clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"`"/Library/Frameworks/R.framework/Resources/bin/Rscript" --vanilla -e "cat(system.file('include', 'src', package = 'StanHeaders'))"`" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/publipha/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/phma.cc -o stan_files/phma.o
...
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (publipha)


```
# ssMousetrack

<details>

* Version: 1.1.5
* GitHub: https://github.com/antcalcagni/ssMousetrack
* Source code: https://github.com/cran/ssMousetrack
* Date/Publication: 2019-01-16 17:00:03 UTC
* Number of recursive dependencies: 54

Run `revdepcheck::revdep_details(, "ssMousetrack")` for more info

</details>

## Newly broken

*   checking whether package ‘ssMousetrack’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/ssMousetrack/new/ssMousetrack.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Installation

### Devel

```
* installing *source* package ‘ssMousetrack’ ...
** package ‘ssMousetrack’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/fit_model_gomp.stan
Wrote C++ file "stan_files/fit_model_gomp.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/fit_model_gomp.cc -o stan_files/fit_model_gomp.o
In file included from stan_files/fit_model_gomp.cc:3:
...
  return ::lgamma_r(x, &sign);
         ~~^
/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/stan/math/prim/fun/lgamma.hpp:85:12: error: no member named 'lgamma_r' in the global namespace
  return ::lgamma_r(x, &sign);
         ~~^
1 warning and 2 errors generated.
make: *** [stan_files/fit_model_gomp.o] Error 1
rm stan_files/fit_model_gomp.cc
ERROR: compilation failed for package ‘ssMousetrack’
* removing ‘/Users/andrew/Downloads/StanHeaders/revdep/checks.noindex/ssMousetrack/new/ssMousetrack.Rcheck/ssMousetrack’


```
### CRAN

```
* installing *source* package ‘ssMousetrack’ ...
** package ‘ssMousetrack’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/fit_model_gomp.stan
Wrote C++ file "stan_files/fit_model_gomp.cc"


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/ssMousetrack/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/fit_model_gomp.cc -o stan_files/fit_model_gomp.o
"/Library/Frameworks/R.framework/Resources/bin/Rscript" -e "source(file.path('..', 'tools', 'make_cc.R')); make_cc(commandArgs(TRUE))" stan_files/fit_model_log.stan
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
* DONE (ssMousetrack)


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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/new/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/visit.cc -o stan_files/visit.o
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


clang++ -arch arm64 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"../inst/include" -I"/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include/src" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DBOOST_MATH_OVERFLOW_ERROR_POLICY=errno_on_error -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/BH/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/Rcpp/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/rstan/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/visit/RcppEigen/include' -I'/Users/andrew/Downloads/StanHeaders/revdep/library.noindex/StanHeaders/old/StanHeaders/include' -I/opt/R/arm64/include   -fPIC  -falign-functions=64 -Wall -g -O2  -Wno-unknown-warning-option -Wno-enum-compare -Wno-ignored-attributes -Wno-unused-local-typedef -Wno-sign-compare -Wno-unneeded-internal-declaration -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable -Wno-infinite-recursion -Wno-unknown-pragmas -Wno-unused-lambda-capture -Wno-deprecated-declarations -Wno-deprecated-builtins -Wno-unused-but-set-variables -ftemplate-backtrace-limit=0 -I/opt/homebrew/opt/libomp/include -c stan_files/visit.cc -o stan_files/visit.o

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
