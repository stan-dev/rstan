\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{R2admb News}
\encoding{UTF-8}

\section{Changes in version 0.7.16 (2017-10-29}{
  \itemize{
    \item handle non-standard report files better; add flag to ignore them
    (i.e. not warn) if requested
    \item internal cleanup
  }
}

\section{Changes in version 0.7.15 (2016-12-06}{
  \itemize{
    \item add \code{covfn} argument to \code{read_params}, warn if
  multiple \code{.std} files found (Koen van Benthem)
    \item small CRAN compatibility fixes
  }
}
\section{Changes in version 0.7.14 (2015-7-20}{
  \itemize{
    \item truncate profile files to 8 characters
    (Edgar Gonzalez)
    \item CRAN tweaks
  }
}

\section{Changes in version 0.7.13 (2015-4-12}{
  \itemize{
    \item CRAN tweaks
  }
}

\section{Changes in version 0.7.12 (2015-4-11}{
  \itemize{
    \item add default location for MacOS ( Fernando Cagua )
    \item various R CMD check/roxygenize updates
    \item CRAN tweaks: package title in Title Case
    \item admb_home goes at beginning of path on Unix
  }
}
\section{Changes in version 0.7.11 (2014-5-01}{
  \itemize{
    \item update vignette, documentation tweaks
  }
}

\section{Changes in version 0.7.10 (2013-9-22}{
  \itemize{
    \item test tweaks
  }
}

\section{Changes in version 0.7.9 (2013-9-21}{
  \itemize{
    \item CRAN tweaks
    \item admb_version now returns NA if ADMB not found
  }
}
  
\section{Changes in version 0.7.8 (2013-9-20}{
  \itemize{
    \item bug fix for AIC
    \item improvements to read_pars
    \item allow all-numeric report file
    \item bug fix for parameter section containing only LOCAL_CALCS
  }
}


\section{Changes in version 0.7.7 (2013-2-20}{
  \itemize{
    \item allow multiple LOC_CALCS blocks per section 
    \item improved handling of variable types (fixed, random, extra) in
  accessors (coef, confint, stdEr)
    \item phase support
    \item add AIC to summary printing method
    \item read_plt() (likelihood profile info) now user-visible
    \item warn on missing cor file
    \item add admb_version()
    \item fix buglets with MCMC reading
    \item fix bug to allow re to be specified within parameter lists
    \item more flexibility in extracting coefficients, var-cov elements
  }
}

\section{Changes in version 0.7.6 (2012-11-06)}{
  \itemize{
    \item Development migration to github
    \item Add Jeff Laake to authors
    \item Improved reading of hessian, vcov matrix
  }
}

\section{Changes in version 0.7.5.3}{
  \itemize{
    \item ignore !! *after* beginning of line in variable definitions
    \item fix mcmc2 bug
  }
}

\section{Changes in version 0.7.5.2}{
  \itemize{
    \item fix Windows bug in run_admb
  }
}

\section{Changes in version 0.7.5.1}{
  \itemize{
    \item export dat_write, add append option, doc tweaks
    \item structure cleanup
  }
}

\section{Changes in version 0.7.5}{
  \itemize{
    \item remove ggplot2 dependency (removed 'ggplot' option from
    plot.admb.hist) [to reduce dependency creep in glmmADMB]
    \item allow for negative-phase in read_par
  }
}

\section{Changes in version 0.7.4}{
  \itemize{
    \item major changes to modularity: exposed compile_admb, run_admb,
    read_admb
    \item changes to clean_admb
    \item improved logic for retrieving profile and MCMC parameter names
    in read_admb; run_admb writes .mcinfo file to disk
    \item write_dat complains if class is not matrix or factor
    \item exposed confint.admb
    \item read_pars reads full std matrix (attempt at
    back-compatibility) [and fixed bug]
    \item fixed bug: std err/estimate names mismatch
    
}}

\section{Changes in version 0.7.2}{
  \itemize{
    \item exposed and documented read_pars, write_pin, write_dat
}}

\section{Changes in version 0.7.1 (2010-08-16)}{
  \itemize{
    \item improved warning handling
    \item took out testing for spurious xxalloc error messages
    \item added ignore_admb_errors
}}

\section{Changes in version 0.7}{
  \itemize{
    \item uncertain: more accessors, general cleanup
}}

\section{Changes in version 0.6 (2010-08-16)}{
  \itemize{
    \item Support for MCMC
    \item Support for profiling
    \item Support for random effects
    \item Expanded vignette
    \item PARAMETERS and DATA section can be created automatically
}}

\section{Changes in version 0.5 (2009-11-25)}{
  \itemize{
    \item Initial release.
    Basic stuff works on my system, lightly tested elsewhere.
}}
