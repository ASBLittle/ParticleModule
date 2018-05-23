node {
   stage('Configure') {
      sh '''#!/bin/bash -l
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .
'''
   }
   stage('Build') {
      sh '''#!/bin/bash -l
make'''
   }
   stage('test') {
      sh '''#!/bin/bash -l
export PYTHONPATH=$PWD
export PATH=$PWD/bin:$PATH
py.test-2.7 --junit-xml=test_results.xml --junit-prefix=Particles \
    --cov --cov-report=html --cov-report=xml
   '''
   }
   stage('Results') {
      junit '**/test_results.xml'
   }
}