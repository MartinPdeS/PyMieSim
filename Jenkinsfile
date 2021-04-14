pipeline {
    agent { docker { image 'python:3.8' } }
    stages {
        stage('build') {
            steps {
                sh 'python --version'
                sh 'pip3 install -r requirements.txt'
                sh 'apt-get install libboost-all-dev'
            }
        }
    }
}
