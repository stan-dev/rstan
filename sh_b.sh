#!/bin/bash

red='\033[0;31m'
NC='\033[0m' # no color

STAN_REPO_BRANCH=develop
grepstanbranch=`git ls-remote --heads https://github.com/stan-dev/stan.git | grep "/${STAN_REPO_BRANCH}"`
if [ -z "$grepstanbranch" ]; then
    echo -e "${red}ERROR:${NC} stan repo does not have {STAN_REPO_BRANCH}"
    exit 20
fi

git config --file=.gitmodules -l
git config -f .gitmodules submodule.stan.branch ${STAN_REPO_BRANCH}
git submodule update --init --recursive --remote --force
git submodule status

rm -Rf StanHeaders/inst/include/upstream StanHeaders/inst/include/src StanHeaders/inst/include/mathlib StanHeaders/inst/include/stan StanHeaders/inst/include/libsundials
cp -Rpv --remove-destination stan/ StanHeaders/inst/include/upstream
cp -Rpv --remove-destination stan/src StanHeaders/inst/include/src
cp -Rpv --remove-destination stan/lib/stan_math StanHeaders/inst/include/mathlib
cp -Rpv --remove-destination stan/lib/stan_math/stan StanHeaders/inst/include/stan
cp -Rpv --remove-destination stan/lib/stan_math/lib/sundials_5.6.1 StanHeaders/inst/include/libsundials

R CMD build StanHeaders/

stanheadtargz=`find StanHeaders*.tar.gz | sort | tail -n 1`

lookforverfile=`tar ztf ${stanheadtargz} | grep stan/version.hpp`

if [ -z "$lookforverfile" ]; then
    echo -e "${red}ERROR:${NC} stan/version.hpp is not found in StanHeaders pkg"
    exit 2
fi

git checkout .gitmodules
# git submodule deinit -f .

R CMD INSTALL ${stanheadtargz}
