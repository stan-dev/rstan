# IgGeneUsage

<details>

* Version: 1.10.0
* GitHub: https://github.com/snaketron/IgGeneUsage
* Source code: https://github.com/cran/IgGeneUsage
* Date/Publication: 2022-04-26
* Number of recursive dependencies: 101

Run `revdep_details(, "IgGeneUsage")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      > test_check("IgGeneUsage")
      Tests Frequentist methods 
      Tests Summarized Experiment check 
      Tests input rules 
      Tests stan model 
      [ FAIL 0 | WARN 0 | SKIP 0 | PASS 91 ]
      > 
      > proc.time()
         user  system elapsed 
       60.855   4.420  65.434 
      
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
* Number of recursive dependencies: 102

Run `revdep_details(, "PosteriorBootstrap")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
       1. ├─base::suppressWarnings(...) at test_anpl.R:177:2
       2. │ └─base::withCallingHandlers(...)
       3. ├─rstan::vb(...)
       4. └─rstan::vb(...)
       5.   └─rstan .local(object, ...)
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
* Number of recursive dependencies: 71

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
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/new/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/BH/include'    -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-kt0bjq/r-base-4.2.0=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -c RcppExports.cpp -o RcppExports.o
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
g++ -std=gnu++14 -I"/usr/share/R/include" -DNDEBUG  -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/Rcpp/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/RcppEigen/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/old/StanHeaders/include' -I'/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/ProbReco/BH/include'    -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-kt0bjq/r-base-4.2.0=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -c RcppExports.cpp -o RcppExports.o
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
* Number of recursive dependencies: 97

Run `revdep_details(, "prophet")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      Chain 4:                0.104 seconds (Sampling)
      Chain 4:                6.651 seconds (Total)
      Chain 4: 
      [ FAIL 0 | WARN 1 | SKIP 0 | PASS 375 ]
      
      [ FAIL 0 | WARN 1 | SKIP 0 | PASS 375 ]
      > 
      > proc.time()
         user  system elapsed 
      153.791   4.536 158.624 
      
       *** caught segfault ***
      address 0x18, cause 'memory not mapped'
      An irrecoverable exception occurred. R is aborting now ...
      Segmentation fault (core dumped)
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 53.7Mb
      sub-directories of 1Mb or more:
        libs  53.0Mb
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

* Version: 2.21.3
* GitHub: https://github.com/stan-dev/rstanarm
* Source code: https://github.com/cran/rstanarm
* Date/Publication: 2022-04-09 00:10:02 UTC
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
      /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rstanarm/RcppEigen/include/Eigen/src/Core/DenseCoeffsBase.h:650:74: warning: ignoring attributes on template argument ‘Eigen::internal::packet_traits<double>::type’ {aka ‘__m128d’} [-Wignored-attributes]
        650 |   return internal::first_aligned<int(unpacket_traits<DefaultPacketType>::alignment),Derived>(m);
            |                                                                          ^~~~~~~~~
      g++ -std=gnu++14 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o sourceCpp_2.so file10bf6559c76c6.o /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rstanarm/rstan/lib//libStanServices.a -L/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/new/StanHeaders/lib/ -lStanHeaders -L/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rstanarm/RcppParallel/lib/ -ltbb -L/usr/lib/R/lib -lR
      [ FAIL 0 | WARN 0 | SKIP 0 | PASS 191 ]
      > 
      > 
      > proc.time()
         user  system elapsed 
       57.007   3.303  60.644 
      
       *** caught segfault ***
      address 0x18, cause 'memory not mapped'
      An irrecoverable exception occurred. R is aborting now ...
      Segmentation fault (core dumped)
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 399.5Mb
      sub-directories of 1Mb or more:
        R       1.6Mb
        doc     4.8Mb
        libs  392.4Mb
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

