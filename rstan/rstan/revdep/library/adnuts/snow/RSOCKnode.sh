#! /bin/sh


# the & for backgrounding works in bash--does it work in other sh variants?

${RPROG:-R} --vanilla <<EOF > ${OUT:-/dev/null} 2>&1 &

library(snow)

options(timeout=getClusterOption("timeout"))

slaveLoop(makeSOCKmaster())
EOF
