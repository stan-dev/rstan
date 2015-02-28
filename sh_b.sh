#!/bin/sh 

STAN_REPO_BRANCH=develop
grepstanbranch=`git ls-remote --heads https://github.com/stan-dev/stan.git | grep "/${STAN_REPO_BRANCH}"`
if [ -z "$grepstanbranch" ]; then
    echo "stan repo does not have {STAN_REPO_BRANCH}"
    exit 20
fi

git config -f .gitmodules submodule.stan.branch ${STAN_REPO_BRANCH}
git submodule update --init --remote
git submodule status

R CMD build StanHeaders/

stanheadtargz=`find StanHeaders*.tar.gz`

lookforverfile=`tar ztf ${stanheadtargz} | grep stan/version.hpp`

if [ -z "$lookforverfile" ]; then
    echo "stan/version.hpp is not found in StanHeaders pkg"
    exit 2
fi

git checkout .gitmodules
# git submodule deinit -f .

R CMD INSTALL ${stanheadtargz}
