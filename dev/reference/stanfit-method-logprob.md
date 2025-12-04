# `log_prob` and `grad_log_prob` functions

Using model's `log_prob` and `grad_log_prob` take values from the
unconstrained space of model parameters and (by default) return values
in the same space. Sometimes we need to convert the values of parameters
from their support defined in the parameters block (which might be
constrained, and for simplicity, we call it the constrained space) to
the unconstrained space and vice versa. The `constrain_pars` and
`unconstrain_pars` functions are used for this purpose.

## Usage

    <!-- %% log_prob(object, upars, adjust_transform = TRUE, gradient = FALSE)  -->
      # S4 method for class 'stanfit'
    log_prob(object, upars, adjust_transform = TRUE, gradient = FALSE)
      <!-- %% grad_log_prob(object, upars, adjust_transform = TRUE)  -->
      # S4 method for class 'stanfit'
    grad_log_prob(object, upars, adjust_transform = TRUE)
      <!-- %% get_num_upars(object) -->
      # S4 method for class 'stanfit'
    get_num_upars(object)
      <!-- %% constrain_pars(object, upars) -->
      # S4 method for class 'stanfit'
    constrain_pars(object, upars)
      <!-- %% unconstrain_pars(object, pars) -->
      # S4 method for class 'stanfit'
    unconstrain_pars(object, pars)

## Methods

- log_prob:

  `signature(object = "stanfit")` Compute `lp__`, the log posterior (up
  to an additive constant) for the model represented by a `stanfit`
  object. Note that, by default, `log_prob` returns the log posterior in
  the *unconstrained* space Stan works in internally. set
  `adjust_transform = FALSE` to make the values match Stan's output.

- grad_log_prob:

  `signature(object = "stanfit")` Compute the gradients for `log_prob`
  as well as the log posterior. The latter is returned as an attribute.

- get_num_upars:

  `signature(object = "stanfit")` Get the number of unconstrained
  parameters.

- constrain_pars:

  `signature(object = "stanfit")` Convert values of the parameter from
  unconstrained space (given as a vector) to their constrained space
  (returned as a named list).

- unconstrain_pars:

  `signature(object = "stanfit")` Contrary to `constrained`, conert
  values of the parameters from constrained to unconstrained space.

## Arguments

- object:

  An object of class
  [`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

- pars:

  An list specifying the values for all parameters on the constrained
  space.

- upars:

  A numeric vector for specifying the values for all parameters on the
  unconstrained space.

- adjust_transform:

  Logical to indicate whether to adjust the log density since Stan
  transforms parameters to unconstrained space if it is in constrained
  space. Set to `FALSE` to make the function return the same values as
  Stan's `lp__` output.

- gradient:

  Logical to indicate whether gradients are also computed as well as the
  log density.

## Details

Stan requires that parameters be defined along with their support. For
example, for a variance parameter, we must define it on the positive
real line. But inside Stan's samplers all parameters defined on the
constrained space are transformed to an unconstrained space amenable to
Hamiltonian Monte Carlo. Because of this, Stan adjusts the log density
function by adding the log absolute value of the Jacobian determinant.
Once a new iteration is drawn, Stan transforms the parameters back to
the original constrained space without requiring interference from the
user. However, when using the log density function for a model exposed
to R, we need to be careful. For example, if we are interested in
finding the mode of parameters on the constrained space, we then do not
need the adjustment. For this reason, the `log_prob` and `grad_log_prob`
functions accept an `adjust_transform` argument.

## Value

`log_prob` returns a value (up to an additive constant) the log
posterior. If `gradient` is `TRUE`, the gradients are also returned as
an attribute with name `gradient`.

`grad_log_prob` returns a vector of the gradients. Additionally, the
vector has an attribute named `log_prob` being the value the same as
`log_prob` is called for the input parameters.

`get_num_upars` returns the number of parameters on the unconstrained
space.

`constrain_pars` returns a list and `unconstrain_pars` returns a vector.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org>.

## See also

[`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md)

## Examples

``` r
# \dontrun{
# see the examples in the help for stanfit as well
# do a simple optimization problem 
opcode <- "
parameters {
  real y;
}
model {
  target += log(square(y - 5) + 1);
}
"
opfit <- stan(model_code = opcode, chains = 0)
#> the number of chains is less than 1; sampling not done
tfun <- function(y) log_prob(opfit, y)
tgrfun <- function(y) grad_log_prob(opfit, y)
or <- optim(1, tfun, tgrfun, method = 'BFGS')
print(or)
#> $par
#> [1] 5
#> 
#> $value
#> [1] 0
#> 
#> $counts
#> function gradient 
#>       30       11 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
#> 

# return the gradient as an attribute
tfun2 <- function(y) { 
  g <- grad_log_prob(opfit, y) 
  lp <- attr(g, "log_prob")
  attr(lp, "gradient") <- g
  lp
} 

or2 <- nlm(tfun2, 10)
or2 
#> $minimum
#> [1] 0
#> 
#> $estimate
#> [1] 5
#> 
#> $gradient
#> [1] 7.897683e-12
#> 
#> $code
#> [1] 1
#> 
#> $iterations
#> [1] 12
#> 
# }
```
