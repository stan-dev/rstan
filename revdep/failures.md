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
