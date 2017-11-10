pipeline {
    agent any
    stages {
        stage('Install dependencies') {
            steps {
                sh """
                    R -e 'install.packages("devtools")'
                    R -e 'update(devtools::package_deps("rstan"))'
                """
            }
        }
        stage('Build') {
            steps {
                sh """
                    export MAKEFLAGS=-j${env.PARALLEL}
                    export CC=${env.CXX}
                    R CMD build rstan/rstan
                """
            }
        }
        stage("Check timings and output") {
            steps {
                sh """
                    R CMD check --as-cran --timings rstan_*.tar.gz || rstan.Rcheck/00install.out
                """
            }
        }
        stage("Check additional unit tests") {
            steps {
                sh """
                    R CMD INSTALL rstan_*.tar.gz
                    cd rstan
                    make test-R
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
