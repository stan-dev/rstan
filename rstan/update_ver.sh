#!/bin/bash
# OLD_VER=[[:digit:]]\\+\\.[[:digit:]]\\+\\.[[:digit:]]

OLD_VER=2.5.0
NEW_VER=2.6.0
SED=sed

${SED} -i "s/${OLD_VER}/${NEW_VER}/g" ./rstan/DESCRIPTION ./rstan/inst/CITATION ./rstan/man/rstan.Rd ./rstan/vignettes/rstan.bib 

OLD_VER2=${OLD_VER//./}
NEW_VER2=${NEW_VER//./}

${SED} -i "s/stanc${OLD_VER2}/stanc${NEW_VER2}/g" ./rstan/R/stanc.R ./rstan/src/init.cpp ./rstan/src/stanc.cpp

echo "use git diff to double check" 
echo "please change dates in DESCRIPTION and man/rstan.md"
