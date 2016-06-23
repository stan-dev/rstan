\name{Plots}
\alias{stan_plot}
\alias{stan_trace}
\alias{stan_dens}
\alias{stan_hist}
\alias{stan_scat}
\alias{stan_ac}
\alias{quietgg}

\title{ggplot2 for RStan}
\description{
  Visual posterior analysis using ggplot2.
}
\usage{
  stan_plot(object, pars, include = TRUE, unconstrain = FALSE, ...)
  stan_trace(object, pars, include = TRUE, unconstrain = FALSE,
            inc_warmup = FALSE, nrow = NULL, ncol = NULL, ..., 
            window = NULL)
  stan_scat(object, pars, unconstrain = FALSE,
            inc_warmup = FALSE, nrow = NULL, ncol = NULL, ...)
  stan_hist(object, pars, include = TRUE, unconstrain = FALSE,
            inc_warmup = FALSE, nrow = NULL, ncol = NULL, ...)
  stan_dens(object, pars, include = TRUE, unconstrain = FALSE,
            inc_warmup = FALSE, nrow = NULL, ncol = NULL, ..., 
            separate_chains = FALSE)
  stan_ac(object, pars, include = TRUE, unconstrain = FALSE,
            inc_warmup = FALSE, nrow = NULL, ncol = NULL, ..., 
            separate_chains = FALSE, lags = 25, partial = FALSE)
  quietgg(gg)
}
\arguments{
  \item{object}{A stanfit or stanreg object.}
  
  \item{pars}{Optional character vector of parameter names. 
    If \code{object} is a stanfit object, the default is to show
    all user-defined parameters or the first 10 (if there are more 
    than 10). If \code{object} is a stanreg object, the default is 
    to show all (or the first 10) regression coefficients 
    (including the intercept). For \code{stan_scat} only, 
    \code{pars} should not be missing and should contain exactly 
    two parameter names.}
  
  \item{include}{Should the parameters given by the \code{pars} argument be
    included (the default) or excluded from the plot?}
  
  \item{unconstrain}{Should parameters be plotted on the unconstrained space?
    Defaults to \code{FALSE}. Only available if \code{object} is a 
    stanfit object.}
  
  \item{inc_warmup}{Should warmup iterations be included? Defaults to
    \code{FALSE}.}

  \item{nrow,ncol}{Passed to \code{\link[ggplot2]{facet_wrap}}.}
  
  \item{...}{Optional additional named arguments passed to geoms 
  (e.g. for \code{stan_trace} the geom is \code{\link[ggplot2]{geom_path}} 
  and we could specify \code{linetype}, \code{size}, \code{alpha}, etc.). 
  For \code{stan_plot} there are also additional arguments that can be specified 
  in \code{...} (see Details).}
  
  \item{window}{For \code{stan_trace} \code{window} is used to control
  which iterations are shown in the plot. See \code{\link{traceplot}}.}
  
  \item{separate_chains}{For \code{stan_dens}, should the density for each 
  chain be plotted? The default is \code{FALSE}, which means that for each 
  parameter the draws from all chains are combined. For \code{stan_ac}, 
  if \code{separate_chains=FALSE} (the default), the autocorrelation is
  averaged over the chains. If \code{TRUE} each chain is plotted separately.}
  
  \item{lags}{For \code{stan_ac}, the maximum number of lags to show.}
  
  \item{partial}{For \code{stan_ac}, should partial autocorrelations be 
  plotted instead? Defaults to \code{FALSE}.}
  
  \item{gg}{A ggplot object or an expression that creates one.}
}
\value{
  A \code{\link[ggplot2]{ggplot}} object that can be further customized
  using the \pkg{ggplot2} package.
}
\note{
Because the \pkg{rstan} plotting functions use \pkg{ggplot2} (and thus the 
resulting plots behave like \code{ggplot} objects), when calling a plotting 
function within a loop or when assigning a plot to a name 
(e.g., \code{graph <- plot(fit, plotfun = "rhat")}),
if you also want the side effect of the plot being displayed you 
must explicity print it (e.g., \code{(graph <- plot(fit, plotfun = "rhat"))}, 
\code{print(graph <- plot(fit, plotfun = "rhat"))}).
}
\details{
For \code{stan_plot}, there are additional arguments that can be specified in 
\code{...}. The optional arguments and their default values are:
\describe{
  \item{\code{point_est = "median"}}{The point estimate to show. 
    Either "median" or "mean".
    }
  \item{\code{show_density = FALSE}}{Should kernel density estimates be plotted 
    above the intervals?
    }
  \item{\code{ci_level = 0.8}}{The posterior uncertainty interval to highlight.
    Central \code{100*ci_level}\% intervals are computed from the quantiles of 
    the posterior draws.}
  \item{\code{outer_level = 0.95}}{An outer interval to also draw as a line 
    (if \code{show_outer_line} is \code{TRUE}) but not highlight.
    }
  \item{\code{show_outer_line = TRUE}}{Should the \code{outer_level} interval 
    be shown or hidden? Defaults to = \code{TRUE} (to plot it).
    }
  \item{\code{fill_color}, \code{outline_color}, \code{est_color}}{Colors to 
    override the defaults for the highlighted interval, the outer interval 
    (and density outline), and the point estimate.
    }
  }
}
\seealso{
\code{\link[=rstan-plotting-functions]{List of RStan plotting functions}},
\code{\link[=rstan_gg_options]{Plot options}}
}
\examples{
example("read_stan_csv")
stan_plot(fit)
stan_trace(fit)

\dontrun{
library(gridExtra)
fit <- stan_demo("eight_schools")

stan_plot(fit)
stan_plot(fit, point_est = "mean", show_density = TRUE, fill_color = "maroon")


# histograms
stan_hist(fit)
# suppress ggplot2 messages about default bindwidth
quietgg(stan_hist(fit))
quietgg(h <- stan_hist(fit, pars = "theta", binwidth = 5)) 

# juxtapose histograms of tau and unconstrained tau 
tau <- stan_hist(fit, pars = "tau")
tau_unc <- stan_hist(fit, pars = "tau", unconstrain = TRUE) + 
            xlab("tau unconstrained")
grid.arrange(tau, tau_unc)

# kernel density estimates
stan_dens(fit)
(dens <- stan_dens(fit, fill = "skyblue", ))
dens <- dens + ggtitle("Kernel Density Estimates\\n") + xlab("")
dens

(dens_sep <- stan_dens(fit, separate_chains = TRUE, alpha = 0.3))
dens_sep + scale_fill_manual(values = c("red", "blue", "green", "black"))
(dens_sep_stack <- stan_dens(fit, pars = "theta", alpha = 0.5,
                             separate_chains = TRUE, position = "stack"))

# traceplot
trace <- stan_trace(fit)
trace +
  scale_color_manual(values = c("red", "blue", "green", "black"))
trace +
  scale_color_brewer(type = "div") +
  theme(legend.position = "none")

facet_style <- theme(strip.background = element_rect(fill = "white"),
                     strip.text = element_text(size = 13, color = "black"))
(trace <- trace + facet_style)

# scatterplot
(mu_vs_tau <- stan_scat(fit, pars = c("mu", "tau"), color = "blue", size = 4))
mu_vs_tau + 
  coord_flip() + 
  theme(panel.background = element_rect(fill = "black"))
  
}
}