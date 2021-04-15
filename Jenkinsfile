pipeline {

    agent { docker { image 'python:3.8'
                     args  '-u root:sudo'}
          }

    stages {
        stage('build') {
            steps {
                sh 'python --version'
                sh 'pip3 install vtk'
                sh 'apt-get update && apt-get install -y python3-opencv'
                sh 'pip3 install opencv-python'
                sh 'pip3 install -r requirements.txt'
                sh 'apt-get install -y libboost-all-dev'
                sh 'python3 setup.py build_ext'
                sh 'pip3 install .'
                sh 'pip3 install codecov'

            }
        }
    }
}
