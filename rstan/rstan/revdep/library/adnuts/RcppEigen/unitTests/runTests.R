
pkg <- "RcppEigen"

if ( ! require( "inline", character.only = TRUE, quietly = TRUE ) ){
    stop( "The inline package is required to run RcppEigen unit tests" )
}

if ( compareVersion( packageDescription( "inline" )[["Version"]], "0.3.5" ) < 0 ){
    stop( "RcppEigen unit tests need at least the version 0.3.5 of inline" )
}

if (require("RUnit", quietly = TRUE)) {

    is_local <- function(){
    	if ( exists( "argv", globalenv() ) && "--local" %in% argv ) return(TRUE)
    	if ( "--local" %in% commandArgs(TRUE) ) return(TRUE)
    	FALSE
    }
    if ( is_local() ) path <- getwd()

    library(package=pkg, character.only = TRUE)
    if (!(exists("path") && file.exists(path)))
        path <- system.file("unitTests", package = pkg)

    ## --- Testing ---

    ## Define tests
    testSuite <- defineTestSuite(name=paste(pkg, "unit testing"), dirs = path)

    if (interactive()) {
        cat("Now have RUnit Test Suite 'testSuite' for package '", pkg, "' :\n", sep='')
        str(testSuite)
        cat('', "Consider doing",
            "\t  tests <- runTestSuite(testSuite)", "\nand later",
            "\t  printTextProtocol(tests)", '', sep="\n")
    } else { ## run from shell / Rscript / R CMD Batch / ...
        ## Run
        tests <- runTestSuite(testSuite)

        output <- NULL

        process_args <- function(argv){
            if ( !is.null(argv) && length(argv) > 0 ){
                rx <- "^--output=(.*)$"
                g  <- grep( rx, argv, value = TRUE )
                if ( length(g) ){
                    sub( rx, "\\1", g[1L] )
                }
            }
        }

        # R CMD check uses this
        if ( exists( "RcppEigen.unit.test.output.dir", globalenv() ) ){
            output <- RcppEigen.unit.test.output.dir
        } else {

            ## give a chance to the user to customize where he/she wants
            ## the unit tests results to be stored with the --output= command
            ## line argument
            if ( exists( "argv",  globalenv() ) ){
                ## littler
                output <- process_args(argv)
            } else {
                ## Rscript
                output <- process_args(commandArgs(TRUE))
            }
        }

        if( is.null(output) ) {         # if it did not work, use parent dir
            output <- ".."              # as BDR does not want /tmp to be used
        }

        ## Print results
        output.txt  <- file.path( output, sprintf("%s-unitTests.txt", pkg))
        output.html <- file.path( output, sprintf("%s-unitTests.html", pkg))

        printTextProtocol(tests, fileName=output.txt)
        message( sprintf( "saving txt unit test report to '%s'", output.txt ) )

        ## Print HTML version to a file
        ## printHTMLProtocol has problems on Mac OS X
        if (Sys.info()["sysname"] != "Darwin"){
            message( sprintf( "saving html unit test report to '%s'", output.html ) )
            printHTMLProtocol(tests, fileName=output.html)
        }

        ##  stop() if there are any failures i.e. FALSE to unit test.
        ## This will cause R CMD check to return error and stop
        err <- getErrors(tests)
        if ( (err$nFail + err$nErr) > 0) {
            stop( sprintf( "unit test problems: %d failures, %d errors", err$nFail, err$nErr) )
        } else {
            success <- err$nTestFunc - err$nFail - err$nErr - err$nDeactivated
            cat( sprintf( "%d / %d\n", success, err$nTestFunc ) )
        }
    }
} else {
    cat("R package 'RUnit' cannot be loaded -- no unit tests run\n", "for package", pkg,"\n")
}

