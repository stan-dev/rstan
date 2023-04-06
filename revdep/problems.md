# prophet

<details>

* Version: 1.0
* GitHub: https://github.com/facebook/prophet
* Source code: https://github.com/cran/prophet
* Date/Publication: 2021-03-30 12:10:07 UTC
* Number of recursive dependencies: 101

Run `revdepcheck::revdep_details(, "prophet")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      Chain 4:                3.463 seconds (Sampling)
      Chain 4:                10.454 seconds (Total)
      Chain 4: 
      [ FAIL 0 | WARN 2 | SKIP 0 | PASS 375 ]
      
      [ FAIL 0 | WARN 2 | SKIP 0 | PASS 375 ]
      > 
      > proc.time()
         user  system elapsed 
      152.305   2.580 154.924 
      
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

* Version: 2.21.3
* GitHub: https://github.com/stan-dev/rstanarm
* Source code: https://github.com/cran/rstanarm
* Date/Publication: 2022-04-09 00:10:02 UTC
* Number of recursive dependencies: 136

Run `revdepcheck::revdep_details(, "rstanarm")` for more info

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
      g++ -std=gnu++14 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o sourceCpp_2.so file3c4d6a434894.o /mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rstanarm/rstan/lib//libStanServices.a -L/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/StanHeaders/new/StanHeaders/lib/ -lStanHeaders -L/mnt/local_drive/badr/tmp/revdepcheck/StanHeaders_revdep/StanHeaders/revdep/library/rstanarm/RcppParallel/lib/ -ltbb -L/usr/lib/R/lib -lR
      [ FAIL 0 | WARN 0 | SKIP 0 | PASS 191 ]
      > 
      > 
      > proc.time()
         user  system elapsed 
       56.130   4.353  60.573 
      
       *** caught segfault ***
      address 0x18, cause 'memory not mapped'
      An irrecoverable exception occurred. R is aborting now ...
      Segmentation fault (core dumped)
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 399.3Mb
      sub-directories of 1Mb or more:
        R       1.6Mb
        doc     4.8Mb
        libs  392.2Mb
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

