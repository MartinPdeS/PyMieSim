pipeline {

    agent { docker { image 'python:3.8'
                     args  '-u root:sudo'}
          }

    stages {
        stage('build') {
            steps {
                sh 'python --version'
                sh 'pip3 install vtk'
                sh 'pip3 install -r requirements.txt'
                sh 'apt-get install libboost-all-dev'
            }
        }
    }
}
