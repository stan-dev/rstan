local({
    master <- "localhost"
    port <- "8765"
    snowlib <- Sys.getenv("R_SNOW_LIB")
    outfile <- Sys.getenv("R_SNOW_OUTFILE")

    args <- commandArgs()
    pos <- match("--args", args)
    args <- args[-(1 : pos)]
    for (a in args) {
        pos <- regexpr("=", a)
        name <- substr(a, 1, pos - 1)
        value <- substr(a,pos + 1, nchar(a))
        switch(name,
               MASTER = master <- value,
               PORT = port <- value,
               SNOWLIB = snowlib <- value,
               OUT = outfile <- value,
               RANK = rank <- value,
               TMPWS = tmpWsName <- value)
    }
    ##**** these should be passed as arguments to makeNWSmaster
    Sys.setenv(MASTER = master)
    Sys.setenv(PORT = port)
    Sys.setenv(RANK = rank)
    Sys.setenv(TMPWS = tmpWsName)

    if (! (snowlib %in% .libPaths()))
        .libPaths(c(snowlib, .libPaths()))
    library(methods) ## because Rscript as of R 2.7.0 doesn't load methods
    library(nws)
    library(snow)

    sinkWorkerOutput(outfile)
    master <- makeNWSmaster()
    sendData(master, "ping")
    cat("starting NWS worker\n")
    slaveLoop(master)
})
