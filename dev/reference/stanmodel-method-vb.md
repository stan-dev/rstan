# Run Stan's variational algorithm for approximate posterior sampling

Approximately draw from a posterior distribution using variational
inference.

This is still considered an experimental feature. We recommend calling
[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md) or
[`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md)
for final inferences and only using `vb` to get a rough idea of the
parameter distributions.

## Usage

    <!-- %% vb(object, \dots)   -->
      # S4 method for class 'stanmodel'
    vb(object, data = list(), pars = NA, include = TRUE,
        seed = sample.int(.Machine$integer.max, 1),
        init = 'random', check_data = TRUE,
        sample_file = tempfile(fileext = '.csv'),
        algorithm = c("meanfield", "fullrank"),
        importance_resampling = FALSE, keep_every = 1,
        ...)

## Methods

- vb:

  `signature(object = "stanmodel")` Call Stan's variational Bayes
  methods for the model defined by S4 class `stanmodel` given the data,
  initial values, etc.

## Arguments

- object:

  An object of class
  [`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md).

- data:

  A named `list` or `environment` providing the data for the model or a
  character vector for all the names of objects used as data. See the
  **Passing data to Stan** section in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- pars:

  If not `NA`, then a character vector naming parameters, which are
  included in the output if `include = TRUE` and excluded if
  `include = FALSE`. By default, all parameters are included.

- include:

  Logical scalar defaulting to `TRUE` indicating whether to include or
  exclude the parameters given by the `pars` argument. If `FALSE`, only
  entire multidimensional parameters can be excluded, rather than
  particular elements of them.

- seed:

  The seed for random number generation. The default is generated from 1
  to the maximum integer supported by R on the machine. Even if multiple
  chains are used, only one seed is needed, with other chains having
  seeds derived from that of the first chain to avoid dependent samples.
  When a seed is specified by a number, `as.integer` will be applied to
  it. If `as.integer` produces `NA`, the seed is generated randomly. The
  seed can also be specified as a character string of digits, such as
  `"12345"`, which is converted to integer.

- init:

  Initial values specification. See the detailed documentation for the
  init argument in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- check_data:

  Logical, defaulting to `TRUE`. If `TRUE` the data will be
  preprocessed; otherwise not. See the **Passing data to Stan** section
  in [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- sample_file:

  A character string of file name for specifying where to write samples
  for *all* parameters and other saved quantities. This defaults to a
  temporary file.

- algorithm:

  Either `"meanfield"` (the default) or `"fullrank"`, indicating which
  variational inference algorithm is used. The `"meanfield"` option uses
  a fully factorized Gaussian for the approximation whereas the
  `fullrank` option uses a Gaussian with a full-rank covariance matrix
  for the approximation. Details and additional references are available
  in the Stan manual.

- importance_resampling:

  Logical scalar (defaulting to `FALSE`) indicating whether to do
  importance resampling to adjust the draws at the optimum to be more
  like draws from the posterior distribution

- keep_every:

  Integer scalar (defaulting to 1) indicating the interval by which to
  thin the draws when `imporance_resampling = TRUE`

- ...:

  Other optional parameters:

  - `iter` (positive `integer`), the maximum number of iterations,
    defaulting to 10000.

  - `grad_samples` (positive `integer`), the number of samples for Monte
    Carlo estimate of gradients, defaulting to 1.

  - `elbo_samples` (positive `integer`), the number of samples for Monte
    Carlo estimate of ELBO (objective function), defaulting to 100.
    (ELBO stands for "the evidence lower bound".)

  - `eta` (`double`), positive stepsize weighting parameter for
    variational inference but is ignored if adaptation is engaged, which
    is the case by default.

  - `adapt_engaged` (`logical`), a flag indicating whether to
    automatically adapt the stepsize, defaulting to `TRUE`.

  - `tol_rel_obj` (positive `double`), the convergence tolerance on the
    relative norm of the objective, defaulting to 0.01.

  - `eval_elbo` (positive `integer`), evaluate ELBO every Nth iteration,
    defaulting to 100.

  - `output_samples` (positive `integer`), number of posterior samples
    to draw and save, defaults to 1000.

  - `adapt_iter` (positive `integer`), the maximum number of iterations
    to adapt the stepsize, defaulting to 50. Ignored if
    `adapt_engaged = FALSE`.

  Refer to the manuals for both CmdStan and Stan for more details.

## Value

An object of
[`stanfit-class`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

## See also

[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md)

The manuals of CmdStan and Stan.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org>.

The Stan Development Team *CmdStan Interface User's Guide*.
<https://mc-stan.org>.

## Examples

``` r
# \dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
#> recompiling to avoid crashing R session
f <- vb(m)
#> Chain 1: ------------------------------------------------------------
#> Chain 1: EXPERIMENTAL ALGORITHM:
#> Chain 1:   This procedure has not been thoroughly tested and may be unstable
#> Chain 1:   or buggy. The interface is subject to change.
#> Chain 1: ------------------------------------------------------------
#> Chain 1: 
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Begin eta adaptation.
#> Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
#> Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
#> Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
#> Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
#> Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
#> Chain 1: Success! Found best value [eta = 1] earlier than expected.
#> Chain 1: 
#> Chain 1: Begin stochastic gradient ascent.
#> Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
#> Chain 1:    100           -0.111             1.000            1.000
#> Chain 1:    200            0.083             1.669            2.338
#> Chain 1:    300            0.031             1.679            1.697
#> Chain 1:    400            0.055             1.369            1.697
#> Chain 1:    500           -0.028             1.689            1.697
#> Chain 1:    600           -0.044             1.469            1.697
#> Chain 1:    700            0.022             1.696            1.697
#> Chain 1:    800            0.002             2.683            2.338
#> Chain 1:    900            0.044             2.490            1.697
#> Chain 1:   1000           -0.019             2.572            2.338
#> Chain 1:   1100           -0.053             2.536            2.338   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1200           -0.028             2.389            1.697   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1300           -0.001             5.153            2.970   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1400            0.062             5.211            2.970   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1500           -0.116             5.068            1.538   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1600            0.107             5.239            2.085   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1700            0.036             5.132            1.983   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1800           -0.024             4.424            1.983   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1900           -0.015             4.384            1.983   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2000           -0.102             4.138            1.538   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2100           -0.129             4.096            1.538   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2200           -0.013             4.884            1.983   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2300            0.046             2.079            1.538   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2400            0.143             2.045            1.538   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2500           -0.014             3.035            1.983   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2600           -0.027             2.875            1.285   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2700           -0.001             4.919            1.285   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2800            0.032             4.772            1.035   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2900            0.059             4.761            1.035   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3000            0.103             4.719            1.035   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3100            0.001            13.965            1.285   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3200           -0.114            13.190            1.035   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3300           -0.011            14.024            1.035   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3400            0.050            14.078            1.216   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3500           -0.095            13.087            1.216   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3600            0.088            13.247            1.521   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3700           -0.048            11.288            1.521   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3800           -0.012            11.483            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3900            0.090            11.551            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4000           -0.030            11.908            2.836   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4100            0.074             2.781            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4200           -0.187             2.820            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4300            0.080             2.192            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4400           -0.017             2.627            2.836   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4500           -0.106             2.558            2.836   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4600           -0.005             4.263            2.977   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4700            0.002             4.309            3.297   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4800           -0.029             4.119            3.297   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4900            0.035             4.189            3.297   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5000           -0.060             3.947            1.828   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5100            0.110             3.961            1.828   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5200           -0.088             4.047            2.259   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5300           -0.085             3.715            1.828   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5400            0.112             3.335            1.765   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5500           -0.000            40.169            1.828   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5600           -0.071            38.355            1.765   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5700           -0.081            38.038            1.581   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5800            0.059            38.169            1.765   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5900           -0.032            38.272            1.765   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6000           -0.014            38.246            1.765   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6100           -0.042            38.159            1.765   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6200            0.073            38.091            1.574   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6300           -0.020            38.547            1.765   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6400            0.117            38.488            1.574   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6500           -0.003             5.622            1.574   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6600           -0.077             5.618            1.574   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6700            0.058             5.839            2.338   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6800            0.056             5.604            1.574   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6900           -0.128             5.461            1.435   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7000           -0.204             5.367            1.435   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7100            0.035             5.983            1.574   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7200           -0.078             5.970            1.450   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7300            0.137             5.668            1.450   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7400            0.018             6.214            1.567   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7500           -0.045             2.302            1.450   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7600            0.052             2.393            1.567   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7700           -0.032             2.422            1.567   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7800           -0.072             2.474            1.567   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7900           -0.070             2.334            1.567   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8000           -0.037             2.387            1.567   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8100           -0.019             1.796            1.450   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8200           -0.080             1.728            1.395   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8300            0.020             2.073            1.395   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8400            0.046             1.467            0.934   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8500            0.014             1.549            0.934   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8600           -0.035             1.504            0.934   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8700           -0.014             1.393            0.934   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8800            0.012             1.558            1.408   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8900            0.048             1.631            1.408   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9000            0.096             1.591            1.408   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9100            0.128             1.522            1.408   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9200            0.020             1.990            1.529   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9300           -0.051             1.626            1.408   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9400            0.034             1.820            1.529   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9500           -0.059             1.755            1.529   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9600           -0.164             1.678            1.529   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9700            0.004             6.124            1.574   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9800            0.011             5.969            1.388   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9900            0.014             5.915            1.388   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   10000           -0.134             5.975            1.388   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1: Informational Message: The maximum number of iterations is reached! The algorithm may not have converged.
#> Chain 1: This variational approximation is not guaranteed to be meaningful.
#> Chain 1: 
#> Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
#> Chain 1: COMPLETED.
# }
```
