pipeline {
    agent { dockerfile true }
    stages {
        stage('Build') {
            steps {
                sh """
                    export MAKEFLAGS=-j${env.PARALLEL}
                    cd StanHeaders
                    git submodule update --init --recursive --remote
                    cd ..
                    R CMD build StanHeaders
                    R CMD build --no-build-vignettes rstan/rstan
                """
            }
        }
        stage("Check timings and output") {
            steps {
                sh """
                    export MAKEFLAGS=-j${env.PARALLEL}
                    R CMD check --as-cran --timings StanHeaders_*.tar.gz || \
                      cat StanHeaders.Rcheck/00check.log
                    R CMD INSTALL StanHeaders_*.tar.gz
                    R CMD check --as-cran --timings --run-donttest --run-dontrun rstan_*.tar.gz || \
                      cat rstan.Rcheck/00check.log
                """
            }
        }
        stage("Check additional unit tests") {
            steps {
                sh """
                    export MAKEFLAGS=-j${env.PARALLEL}
                    R CMD INSTALL rstan_*.tar.gz
                    cd rstan
                    make test-R || echo "extra unit tests failed"
                    cd ..
                """
            }
        }
        stage("Check rstanarm") {
            steps {
                sh """
                    export MAKEFLAGS=-j${env.PARALLEL}
                    R -e 'update(devtools::package_deps("rstanarm"), repos = "https://cran.r-project.org")'
                    wget -Nc https://cran.r-project.org/src/contrib/rstanarm_2.17.4.tar.gz
                    R CMD check --as-cran --timings --run-donttest --run-dontrun rstanarm_*.tar.gz || \
                      cat rstanarm.Rcheck/00check.log
                """
            }
        }
    }
    post {
        always {
            warnings consoleParsers: [[parserName: 'GNU C Compiler (gcc)']], failedTotalAll: '0', usePreviousBuildAsReference: false, canRunOnFailed: true
            warnings consoleParsers: [[parserName: 'Clang (LLVM based)']], failedTotalAll: '0', usePreviousBuildAsReference: false, canRunOnFailed: true
            deleteDir()
        }
    }
}
