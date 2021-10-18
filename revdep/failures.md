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
