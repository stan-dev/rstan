#! /bin/sh

${RPROG:-R} --vanilla <<EOF > ${OUT:-/dev/null} 2>&1

library(rpvm)
library(snow)

slaveLoop(makePVMmaster())
.PVM.exit()
EOF
