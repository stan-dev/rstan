#!/bin/bash
# OLD_VER=[[:digit:]]\\+\\.[[:digit:]]\\+\\.[[:digit:]]

OLD_VER=2.2.0
NEW_VER=2.3.0
SED=sed

${SED} -i "s/${OLD_VER}/${NEW_VER}/g" ./rstan/DESCRIPTION ./rstan/inst/CITATION ./rstan/man/rstan.Rd ./rstan/vignettes/rstan.bib 

echo "use git diff to double check" 

