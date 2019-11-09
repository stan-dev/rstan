#! /bin/sh

${RPROG:-R} --vanilla <<EOF > ${OUT:-/dev/null} 2>&1 &

library(nws)
library(snow)

local({
    master <- makeNWSmaster()
    sendData(master, "ping")
    slaveLoop(master)
})
EOF
