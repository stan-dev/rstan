#! /bin/sh

${RPROG:-R} --vanilla <<EOF > ${OUT:-/dev/null} 2>&1

library(Rmpi)
library(snow)

runMPIslave()
EOF
