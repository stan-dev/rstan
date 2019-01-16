FROM r-base

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
  libssl-dev \
  libcurl4-openssl-dev
RUN R -e 'install.packages("devtools", repos = "https://cran.r-project.org")'
RUN R -e 'devtools::install_dev_deps("rstan", repos = "https://cran.r-project.org")'
RUN R -e 'install.packages("RInside", repos = "https://cran.r-project.org")'
