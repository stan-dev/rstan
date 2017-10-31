pipeline {
    agent any
    stages {
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
                    R CMD check --as-cran --timings rstan_*.tar.gz
                    cat rstan.Rcheck/00install.out
                    cat rstan.Rcheck/tests/doRUnit.Rout
                """
            }
        }
    }
    post {
        always {
            warnings consoleParsers: [[parserName: 'GNU C Compiler 4 (gcc)']], failedTotalAll: '0', usePreviousBuildAsReference: false, canRunOnFailed: true
            warnings consoleParsers: [[parserName: 'Clang (LLVM based)']], failedTotalAll: '0', usePreviousBuildAsReference: false, canRunOnFailed: true
            deleteDir()
        }
    }
}
