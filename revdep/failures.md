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
    See ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/nlmixr2est/new/nlmixr2est.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘nlmixr2est’ ...
** package ‘nlmixr2est’ successfully unpacked and MD5 sums checked
** using staged installation
--------[begin src/Makevars]--------
# -*- mode: makefile-gmake -*-
ARMA=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/RcppArmadillo/include
BH=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/BH/include
RCPP=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/Rcpp/include
EG=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/RcppEigen/include
SH=-isystem'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/new/StanHeaders/include/src' 
...
** R
** inst
** byte-compile and prepare package for lazy loading
Error: .onLoad failed in loadNamespace() for 'rxode2parse', details:
  call: NULL
  error: rxode2parse compiled with dparser '1.3.1.9' but dparser '1.3.1.10' is loaded
Recompile rxode2parse with the this version of dparser
Execution halted
ERROR: lazy loading failed for package ‘nlmixr2est’
* removing ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/nlmixr2est/new/nlmixr2est.Rcheck/nlmixr2est’


```
### CRAN

```
* installing *source* package ‘nlmixr2est’ ...
** package ‘nlmixr2est’ successfully unpacked and MD5 sums checked
** using staged installation
--------[begin src/Makevars]--------
# -*- mode: makefile-gmake -*-
ARMA=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/RcppArmadillo/include
BH=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/BH/include
RCPP=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/Rcpp/include
EG=/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/nlmixr2est/RcppEigen/include
SH=-isystem'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/old/StanHeaders/include/src' 
...
** R
** inst
** byte-compile and prepare package for lazy loading
Error: .onLoad failed in loadNamespace() for 'rxode2parse', details:
  call: NULL
  error: rxode2parse compiled with dparser '1.3.1.9' but dparser '1.3.1.10' is loaded
Recompile rxode2parse with the this version of dparser
Execution halted
ERROR: lazy loading failed for package ‘nlmixr2est’
* removing ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/nlmixr2est/old/nlmixr2est.Rcheck/nlmixr2est’


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
    See ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/ProbReco/new/ProbReco.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking installed package size ... NOTE
    ```
      installed size is 19.0Mb
      sub-directories of 1Mb or more:
        libs  18.6Mb
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
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/new/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/BH/include'    -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-ZLat0n/r-base-4.2.2.20221110=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -c RcppExports.cpp -o RcppExports.o
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
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/old/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/BH/include'    -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-ZLat0n/r-base-4.2.2.20221110=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -c RcppExports.cpp -o RcppExports.o
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
# rxode2ll

<details>

* Version: 2.0.9
* GitHub: https://github.com/nlmixr2/rxode2ll
* Source code: https://github.com/cran/rxode2ll
* Date/Publication: 2022-09-29 08:20:02 UTC
* Number of recursive dependencies: 50

Run `revdepcheck::revdep_details(, "rxode2ll")` for more info

</details>

## Newly broken

*   checking whether package ‘rxode2ll’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/rxode2ll/new/rxode2ll.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking installed package size ... NOTE
    ```
      installed size is 72.9Mb
      sub-directories of 1Mb or more:
        libs  72.7Mb
    ```

## Installation

### Devel

```
* installing *source* package ‘rxode2ll’ ...
** package ‘rxode2ll’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/new/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/BH/include'   -DBOOST_DISABLE_ASSERTS -DBOOST_NO_CXX11_STATIC_ASSERT -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION -DBOOST_NO_AUTO_PTR -D_REENTRANT -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-ZLat0n/r-base-4.2.2.20221110=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -c RcppExports.cpp -o RcppExports.o
In file included from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/Eigen/Core:397,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/Eigen/Dense:1,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/RcppEigenForward.h:30,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/RcppEigen.h:25,
                 from RcppExports.cpp:4:
...
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
Error: package or namespace load failed for ‘rxode2ll’ in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/rxode2ll/new/rxode2ll.Rcheck/00LOCK-rxode2ll/00new/rxode2ll/libs/rxode2ll.so':
  /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/rxode2ll/new/rxode2ll.Rcheck/00LOCK-rxode2ll/00new/rxode2ll/libs/rxode2ll.so: undefined symbol: _ZN3tbb8internal26task_scheduler_observer_v37observeEb
Error: loading failed
Execution halted
ERROR: loading failed
* removing ‘/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/checks/rxode2ll/new/rxode2ll.Rcheck/rxode2ll’


```
### CRAN

```
* installing *source* package ‘rxode2ll’ ...
** package ‘rxode2ll’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/old/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/BH/include'   -DBOOST_DISABLE_ASSERTS -DBOOST_NO_CXX11_STATIC_ASSERT -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION -DBOOST_NO_AUTO_PTR -D_REENTRANT -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-ZLat0n/r-base-4.2.2.20221110=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -c RcppExports.cpp -o RcppExports.o
In file included from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/Eigen/Core:397,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/Eigen/Dense:1,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/RcppEigenForward.h:30,
                 from /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rxode2ll/RcppEigen/include/RcppEigen.h:25,
                 from RcppExports.cpp:4:
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
* DONE (rxode2ll)


```
