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
#> Chain 1: Gradient evaluation took 3e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 1:    100           -0.087             1.000            1.000
#> Chain 1:    200            0.115             1.380            1.759
#> Chain 1:    300            0.011             4.155            1.759
#> Chain 1:    400            0.052             3.314            1.759
#> Chain 1:    500           -0.020             3.377            1.759
#> Chain 1:    600           -0.027             2.861            1.759
#> Chain 1:    700            0.112             2.630            1.244
#> Chain 1:    800           -0.128             2.536            1.759
#> Chain 1:    900           -0.187             2.289            1.244
#> Chain 1:   1000            0.039             2.637            1.759
#> Chain 1:   1100            0.062             2.573            1.759   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1200            0.003             4.084            1.875   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1300            0.045             3.205            1.244   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1400           -0.060             3.302            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1500            0.049             3.162            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1600            0.052             3.141            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1700            0.087             3.056            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1800           -0.019             3.423            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1900           -0.103             3.473            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2000           -0.025             3.203            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2100            0.049             3.318            1.756   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2200            0.108             1.687            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2300            0.072             1.645            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2400           -0.038             1.757            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2500            0.008             2.137            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2600           -0.058             2.244            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2700            0.129             2.349            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2800           -0.015             2.749            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2900            0.141             2.778            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3000            0.018             3.153            1.516   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3100           -0.058             3.132            1.451   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3200            0.188             3.208            1.451   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3300           -0.025             3.997            2.871   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3400            0.038             3.877            1.673   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3500           -0.021             3.554            1.673   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3600           -0.002             4.636            2.802   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3700           -0.117             4.589            2.802   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3800           -0.096             3.656            1.673   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3900           -0.059             3.608            1.673   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4000            0.126             3.074            1.471   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4100           -0.137             3.135            1.673   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4200           -0.153             3.015            1.673   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4300           -0.126             2.197            1.471   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4400            0.007             3.908            1.471   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4500           -0.192             3.732            1.037   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4600            0.038             3.143            1.037   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4700            0.055             3.075            1.037   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4800            0.064             3.068            1.037   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4900            0.039             3.071            1.037   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5000           -0.102             3.062            1.037   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5100            0.135             3.046            1.037   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5200           -0.036             3.507            1.382   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5300           -0.008             3.835            1.752   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5400            0.130             2.062            1.382   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5500           -0.057             2.287            1.752   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5600            0.061             1.874            1.752   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5700           -0.135             1.988            1.752   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5800            0.083             2.237            1.935   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5900           -0.040             2.477            2.631   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6000            0.029             2.579            2.631   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6100           -0.063             2.549            2.631   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6200           -0.003             3.957            2.631   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6300           -0.161             3.704            2.401   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6400           -0.043             3.875            2.631   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6500           -0.167             3.621            2.401   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6600            0.110             3.680            2.522   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6700            0.013             4.291            2.631   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6800            0.026             4.078            2.522   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6900           -0.022             3.988            2.401   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7000            0.098             3.870            2.151   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7100           -0.045             4.045            2.522   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7200            0.060             2.341            2.151   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7300            0.004             3.585            2.522   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7400            0.022             3.390            2.151   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7500           -0.094             3.439            2.151   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7600           -0.035             3.356            1.747   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7700           -0.190             2.682            1.699   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7800           -0.096             2.730            1.699   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7900            0.047             2.819            1.699   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8000           -0.038             2.918            1.747   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8100            0.131             2.728            1.699   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8200            0.029             2.901            1.699   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8300            0.018             1.619            1.293   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8400            0.078             1.615            1.293   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8500           -0.029             1.855            1.699   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8600            0.128             1.808            1.293   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8700           -0.110             1.943            2.168   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8800           -0.019             2.321            2.222   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8900            0.095             2.137            2.168   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9000            0.019             2.323            2.168   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9100           -0.038             2.343            2.168   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9200           -0.059             2.030            1.496   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9300            0.041             2.214            2.168   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9400            0.050             2.155            2.168   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9500            0.096             1.840            1.496   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9600           -0.018             2.339            2.168   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9700           -0.092             2.202            1.496   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9800            0.093             1.925            1.496   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9900           -0.007             3.180            1.984   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   10000            0.168             2.877            1.496   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1: Informational Message: The maximum number of iterations is reached! The algorithm may not have converged.
#> Chain 1: This variational approximation is not guaranteed to be meaningful.
#> Chain 1: 
#> Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
#> Chain 1: COMPLETED.
# }
```
