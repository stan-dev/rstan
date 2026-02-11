# Interfacing with External C++ Code

Starting with the 2.13 release, it is much easier to use external C++
code in a Stan program. This vignette briefly illustrates how to do so.

Suppose that you have (part of) a Stan program that involves Fibonacci
numbers, such as

``` stan
functions {
  int fib(int n);
  int fib(int n) {
    if (n <= 0) reject("n must be positive");
    return n <= 2 ? 1 : fib(n - 1) + fib(n - 2);
  }
}
model {} // use the fib() function somehow
```

On the second line, we have *declared* the `fib` function before it is
*defined* in order to call it recursively.

For functions that are not recursive, it is not necessary to declare
them before defining them but it may be advantageous. For example, I
often like to hide the definitions of complicated utility functions that
are just a distraction using the `#include "file"` mechanism

``` stan
functions {
  real complicated(real a, real b, real c, real d, real e, real f, real g);
#include "complicated.stan"
}
model {} // use the complicated() function somehow
```

This Stan program would have to be parsed using the `stanc_builder`
function in the **rstan** package rather than the default `stanc`
function (which is called by `sampling` and `stan` internally).

Returning to the Fibonacci example, it is not necessary to define the
`fib` function using the Stan language because Stan programs with
functions that are *declared* but not *defined* can use the standard
capabilities of the C++ toolchain to provide the function definitions in
C++. For example, this program produces a parser error by default

``` r
mc <- 
'
functions { int fib(int n); }
model {} // use the fib() function somehow
'
try(stan_model(model_code = mc, model_name = "parser_error"), silent = TRUE)
```

However, if we specify the `allow_undefined` and `includes` arguments to
the `stan_model` function, and define a `fib` function in the named C++
header file, then it will parse and compile

``` r
stan_model(model_code = mc, model_name = "external", allow_undefined = TRUE,
           includes = paste0('\n#include "', 
                             file.path(getwd(), 'fib.hpp'), '"\n'))
```

Specifying the `includes` argument is a bit awkward because the C++
representation of a Stan program is written and compiled in a temporary
directory. Thus, the `includes` argument must specify a *full* path to
the fib.hpp file, which in this case is in the working directory. Also,
the path must be enclosed in double-quotes, which is why single quotes
are used in the separate arguments to the `paste0` function so that
double-quotes are interpreted literally. Finally, the `includes`
argument should include newline characters (`"\n"`) at the start and
end. It is possible to specify multiple paths using additional newline
characters or include a “meta-header” file that contains `#include`
directives to other C++ header files.

The result of the `includes` argument is inserted into the C++ file
directly at the end of the lines (as opposed to CmdStan where it is
inserted directly *before* the following lines)

``` cpp
#include <stan/model/model_header.hpp>

namespace some_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

static int current_statement_begin__;

// various function declarations and / or definitions
#include "/full/path/to/fib.hpp"
```

Thus, the definition of the `fib` function in the fib.hpp file need not
be enclosed in any particular namespace (which is a random string by
default. The “meta-include” stan/model/model_header.hpp file reads as

    #ifndef STAN_MODEL_MODEL_HEADER_HPP
    #define STAN_MODEL_MODEL_HEADER_HPP

    #include <stan/model/model_base.hpp>
    #include <stan/model/model_base_crtp.hpp>
    #include <stan/math.hpp>

    #include <stan/io/deserializer.hpp>
    #include <stan/io/serializer.hpp>

    #include <stan/model/rethrow_located.hpp>
    #include <stan/model/prob_grad.hpp>
    #include <stan/model/indexing.hpp>
    #include <stan/services/util/create_rng.hpp>

    #include <cmath>
    #include <cstddef>
    #include <fstream>
    #include <iostream>
    #include <sstream>
    #include <stdexcept>
    #include <utility>
    #include <vector>

    #endif

so the definition of the `fib` function in the fib.hpp file could
utilize any function in the Stan Math Library (without having to prefix
function calls with `stan::math::`), some typedefs to classes in the
Eigen matrix algebra library, plus streams, exceptions, etc. without
having to worry about the corresponding header files. Nevertheless, an
external C++ file *may* contain additional include directives that bring
in class definitions, for example.

Now let’s examine the fib.hpp file, which contains the C++ lines

``` cpp
int fib(const int&n, std::ostream* pstream__) {
  if (n <= 0) {
    stringstream errmsg;
    errmsg << "n must be positive";
    throw std::domain_error(errmsg.str());
  }
  return n <= 1 ? 1 : fib(n - 1, 0) + fib(n - 2, 0);
}
```

This C++ function is essentially what the preceding user-defined
function in the Stan language

``` stan
int fib(int n) {
  if (n <= 0) reject("n must be positive");
  return n <= 2 ? 1 : fib(n - 1) + fib(n - 2);
}
```

parses to. Thus, there is no *speed* advantage to defining the `fib`
function in the external fib.hpp file. However, it is possible to use an
external C++ file to handle the gradient of a function analytically as
opposed to using Stan’s autodifferentiation capabilities, which are
slower and more memory intensive but fully general. In this case, the
`fib` function only deals with integers so there is nothing to take the
derivative of. The primary advantage of using an external C++ file is
flexibility to do things that cannot be done directly in the Stan
language. It is also useful for R packages like **rstanarm** that may
want to define some C++ functions in the package’s src directory and
rely on the linker to make them available in its Stan programs, which
are compiled at (or before) installation time.

In the C++ version, we check if `n` is non-positive, in which case we
throw an exception. The way the MCMC samplers are configured, if there
is an exception thrown and it is `std::domain_error`, it will treat that
iteration as a rejection but as a recoverable error: set that
iteration’s log probability value to negative infinity and move to the
next iteration. If there is a `std::invalid_argument` thrown, MCMC
terminates. We treat these as non-recoverable, user programming errors.
It is unnecessary to prefix `stringstream` with `std::` because of the
`using std::stringstream;` line in the *generated* C++ file. However,
there is no corresponding `using std::domain_error;` line, so it has to
be qualified appropriately when the exception is thrown.

The only confusing part of the C++ version of the `fib` function is that
it has an additional argument (with no default value) named `pstream__`
that is added to the *declaration* of the `fib` function by Stan. Thus,
your *definition* of the `fib` function needs to match with this
signature. This additional argument is a pointer to a `std::ostream` and
is only used if your function prints something to the screen, which is
rare. Thus, when we call the `fib` function recursively in the last
line, we specify `fib(n - 1, 0) + fib(n - 2, 0);` so that the output (if
any, and in this case there is none) is directed to the null pointer.

This vignette has employed a toy example with the Fibonacci function,
which has little apparent use in a Stan program and if it were useful,
would more easily be implemented as a user-defined function in the
`functions` block as illustrated at the outset. The ability to use
external C++ code only becomes useful with more complicated C++
functions. It goes without saying that this mechanism ordinarily cannot
call functions in C, Fortran, R, or other languages because Stan needs
the derivatives with respect to unknown parameters in order to perform
estimation. These derivatives are handled with custom C++ types that
cannot be processed by functions in other languages that only handle
primitive types such as `double`, `float`, etc.

That said, it is possible to accomplish a great deal in C++,
particularly when utilizing the Stan Math Library. For more details, see
[The Stan Math Library: Reverse-Mode Automatic Differentiation in
C++](https://arxiv.org/abs/1509.07164) and its GitHub
[repository](https://github.com/stan-dev/math/). The functions that you
*declare* in the `functions` block of a Stan program will typically
involve templating and type promotion in their signatures when parsed to
C++ (the only exceptions are functions whose only arguments are
integers, as in the `fib` function above). Suppose you wanted to define
a function whose arguments are real numbers (or at least one of the
arguments is). For example,

``` r
mc <- 
'
functions { real sinc(real x); }
transformed data { real sinc_pi = sinc(pi()); }
'
stan_model(model_code = mc, model_name = "external", allow_undefined = TRUE,
           includes = paste0('\n#include "', 
                             file.path(getwd(), 'sinc.hpp'), '"\n'))
```

    ## Trying to compile a simple C file

    ## Running /opt/R/4.5.2/lib/R/bin/R CMD SHLIB foo.c
    ## using C compiler: ‘gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0’
    ## gcc -std=gnu2x -I"/opt/R/4.5.2/lib/R/include" -DNDEBUG   -I"/home/runner/work/_temp/Library/Rcpp/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/unsupported"  -I"/home/runner/work/_temp/Library/BH/include" -I"/home/runner/work/_temp/Library/StanHeaders/include/src/"  -I"/home/runner/work/_temp/Library/StanHeaders/include/"  -I"/home/runner/work/_temp/Library/RcppParallel/include/"  -I"/home/runner/work/_temp/Library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1  -I/usr/local/include    -fpic  -g -O2  -c foo.c -o foo.o
    ## In file included from <command-line>:
    ## /home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:3:10: fatal error: stdexcept: No such file or directory
    ##     3 | #include <stdexcept>
    ##       |          ^~~~~~~~~~~
    ## compilation terminated.
    ## make: *** [/opt/R/4.5.2/lib/R/etc/Makeconf:202: foo.o] Error 1

    ## Error in `compileCode()`:
    ## !  int Option = 0; Scalar = double]’
    ## /home/runner/work/_temp/Library/StanHeaders/include/src/stan/mcmc/hmc/hamiltonians/dense_e_metric.hpp:22:0:   required from ‘double stan::mcmc::dense_e_metric<Model, BaseRNG>::T(stan::mcmc::dense_e_point&) [with Model = model27d42d94d146_external_namespace::model27d42d94d146_external; BaseRNG = boost::random::mixmax_engine<17, 36, 0>]’
    ## /home/runner/work/_temp/Library/StanHeaders/include/src/stan/mcmc/hmc/hamiltonians/dense_e_metric.hpp:21:0:   required from here
    ## /home/runner/work/_temp/Library/RcppEigen/include/Eigen/src/Core/DenseCoeffsBase.h:654:74: warning: ignoring attributes on template argument ‘Eigen::internal::packet_traits<double>::type’ {aka ‘__m128d’} [-Wignored-attributes]
    ##   654 |   return internal::first_aligned<int(unpacket_traits<DefaultPacketType>::alignment),Derived>(m);
    ##       |                                                                          ^~~~~~~~~
    ## make: *** [/opt/R/4.5.2/lib/R/etc/Makeconf:211: file27d41201ade0.o] Error 1

    ## Error in `sink()`:
    ## ! invalid connection

The sinc.hpp file reads as

    double
    sinc(const double& x, std::ostream* pstream__) {
      return x != 0.0 ? sin(x) / x : 1.0;
    }

    var
    sinc(const var& x, std::ostream* pstream__) {
      double x_ = x.val();
      double f = x_ != 0.0 ? sin(x_) / x_ : 1.0;
      double dfdx_ = x_ != 0.0 ? (cos(x_) - sin(x_)) / x_ : 0.0;
      return var(new precomp_v_vari(f, x.vi_, dfdx_));
    }

The body of the first `sinc` function is simply its mathematical
definition in the form of a ternary operator, which is sufficient when
the input is a `double`.

The last lines of sinc.hpp are a specialization for when the input is an
unknown real parameter, which is represented in the Stan Math Library as
a `stan::math::var`. Since the derivative of the `sinc` function is easy
to compute analytically, we extract the underlying double-precision
value from the inputted `stan::math::var` and use that to calculate the
function value and its first derivative. Then, we return the result of
`precomputed_gradients`, whose arguments are the function value (`f`),
the derivative of `x` with respect to any other parameters (`x.vi_`),
and the first derivative of `f` (`dfdx_`). The latter two are actually
`std::vector`s but only have one element each because there is only one
unknown.

An easy way to see what the generated function signature will be is to
call the `stanc` function in the **rstan** package with
`allow_undefined = TRUE` and inspect the resuling C++ code. In this
case, I first did

``` r
try(readLines(stanc(model_code = mc, allow_undefined = TRUE)$cppcode))
```

    ## Warning in file(con, "r"): expanded path length 13646 would be too long for
    ## #ifndef USE_STANC3
    ## #define USE_STANC3
    ## #endif
    ## // Code generated by stanc v2.38.0-34-g1df15b1
    ## #include <stan/model/model_header.hpp>
    ## namespace model27d451fafaf1_file27d45280c27b_namespace {
    ## using stan::model::model_base_crtp;
    ## using namespace stan::math;
    ## stan::math::profile_map profiles__;
    ## static constexpr std::array<const char*, 2> locations_array__ =
    ##   {" (found before start of program)",
    ##   " (in 'file27d45280c27b', line 2, column 19 to column 45)"};
    ## class model27d451fafaf1_file27d45280c27b final : public model_base_crtp<model27d451fafaf1_file27d45280c27b> {
    ## private:
    ##   double sinc_pi;
    ## public:
    ##   ~model27d451fafaf1_file27d45280c27b() {}
    ##   model27d451fafaf1_file27d45280c27b(stan::io::var_context& context__,
    ##                                      unsigned int random_seed__ = 0,
    ##                                      std::ostream* pstream__ = nullptr)
    ##       : model_base_crtp(0) {
    ##     int current_statement__ = 0;
    ##     // suppress unused var warning [... truncated]
    ## Warning in file(con, "r"): expanded path length 13646 would be too long for
    ## #ifndef USE_STANC3
    ## #define USE_STANC3
    ## #endif
    ## // Code generated by stanc v2.38.0-34-g1df15b1
    ## #include <stan/model/model_header.hpp>
    ## namespace model27d451fafaf1_file27d45280c27b_namespace {
    ## using stan::model::model_base_crtp;
    ## using namespace stan::math;
    ## stan::math::profile_map profiles__;
    ## static constexpr std::array<const char*, 2> locations_array__ =
    ##   {" (found before start of program)",
    ##   " (in 'file27d45280c27b', line 2, column 19 to column 45)"};
    ## class model27d451fafaf1_file27d45280c27b final : public model_base_crtp<model27d451fafaf1_file27d45280c27b> {
    ## private:
    ##   double sinc_pi;
    ## public:
    ##   ~model27d451fafaf1_file27d45280c27b() {}
    ##   model27d451fafaf1_file27d45280c27b(stan::io::var_context& context__,
    ##                                      unsigned int random_seed__ = 0,
    ##                                      std::ostream* pstream__ = nullptr)
    ##       : model_base_crtp(0) {
    ##     int current_statement__ = 0;
    ##     // suppress unused var warning [... truncated]

    ## Warning in file(con, "r"): cannot open file '#ifndef USE_STANC3
    ## #define USE_STANC3
    ## #endif
    ## // Code generated by stanc v2.38.0-34-g1df15b1
    ## #include <stan/model/model_header.hpp>
    ## namespace model27d451fafaf1_file27d45280c27b_namespace {
    ## using stan::model::model_base_crtp;
    ## using namespace stan::math;
    ## stan::math::profile_map profiles__;
    ## static constexpr std::array<const char*, 2> locations_array__ =
    ##   {" (found before start of program)",
    ##   " (in 'file27d45280c27b', line 2, column 19 to column 45)"};
    ## class model27d451fafaf1_file27d45280c27b final : public model_base_crtp<model27d451fafaf1_file27d45280c27b> {
    ## private:
    ##   double sinc_pi;
    ## public:
    ##   ~model27d451fafaf1_file27d45280c27b() {}
    ##   model27d451fafaf1_file27d45280c27b(stan::io::var_context& context__,
    ##                                      unsigned int random_seed__ = 0,
    ##                                      std::ostream* pstream__ = nullptr)
    ##       : model_base_crtp(0) {
    ##     int current_statement__ = 0;
    ##     // suppress unused var warning
    ##     (void) current_statement__; [... truncated]

    ## Error in file(con, "r") : cannot open the connection

to see what function signatures needed to be written for sinc.hpp.

Once you go to the trouble of writing a numerically stable C++ function,
we would welcome a pull request on GitHub to include your C++ function
in the Stan Math Library for everyone to benefit from, provided that it
can be licensed under the 3-clause BSD license and its use is not
overly-specific to your particular Stan program.

The Stan Math Library is compliant with the C++14 standard but not all
compilers fully support the C++14 standard. In particular, the compiler
that comes with RTools does not support all C++14 features. But you can
use many C++14 features, such as the `auto` keyword.
