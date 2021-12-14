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

rm -Rf StanHeaders/inst/include/src \
       StanHeaders/inst/include/mathlib \
       StanHeaders/inst/include/stan \
       StanHeaders/inst/include/libsundials || true

cp -Rf StanHeaders/inst/include/upstream/src StanHeaders/inst/include/src || true
cp -Rf StanHeaders/inst/include/upstream/lib/stan_math StanHeaders/inst/include/mathlib || true
cp -Rf StanHeaders/inst/include/upstream/lib/stan_math/stan StanHeaders/inst/include/stan || true
cp -Rf StanHeaders/inst/include/upstream/lib/stan_math/lib/sundials_* StanHeaders/inst/include/libsundials || true
cp -Rf StanHeaders/inst/include/upstream/lib/stan_math/lib/sundials_*/include/* StanHeaders/inst/include/ || true
cp -Rf StanHeaders/inst/include/upstream/lib/stan_math/lib/sundials_*/src/* StanHeaders/src/ || true

# Prepare TBB objects and headers
if [ -z "$TBB_LIB" ] && [ ! -d StanHeaders/src/tbb ]; then
    cp -Rv StanHeaders/inst/include/upstream/lib/stan_math/lib/tbb_*/src/tbb* StanHeaders/src/ || true
    cp -Rv StanHeaders/inst/include/upstream/lib/stan_math/lib/tbb_*/src/rml StanHeaders/src/ || true
    cp -Rv StanHeaders/inst/include/upstream/lib/stan_math/lib/tbb_*/include/serial StanHeaders/inst/include/ || true
    cp -Rv StanHeaders/inst/include/upstream/lib/stan_math/lib/tbb_*/include/*tbb* StanHeaders/inst/include/ || true
    cp -v StanHeaders/inst/include/upstream/lib/stan_math/lib/tbb_*/include/*tbb*/*.h StanHeaders/inst/include/tbb/ || true
    cp -v StanHeaders/inst/include/upstream/lib/stan_math/lib/tbb_*/src/tbb/*.h StanHeaders/inst/include/tbb/ || true
    echo '#define __TBB_VERSION_STRINGS(N) "Empty"' > StanHeaders/src/tbb/version_string.ver || true
fi

R CMD build "$@" StanHeaders/

stanheadtargz=`find StanHeaders*.tar.gz | sort | tail -n 1`

lookforverfile=`tar ztf ${stanheadtargz} | grep stan/version.hpp`

if [ -z "$lookforverfile" ]; then
    echo -e "${red}ERROR:${NC} stan/version.hpp is not found in StanHeaders pkg"
    exit 2
fi

git checkout .gitmodules
# git submodule deinit -f .

R CMD INSTALL ${stanheadtargz} || Rscript -e 'remotes::install_local(rev(list.files(pattern = Sys.glob("StanHeaders")))[1], dependencies = TRUE, type = "source")'
