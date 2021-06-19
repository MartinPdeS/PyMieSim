FROM quay.io/pypa/manylinux2014_x86_64

ENV QT_QPA_PLATFORM minimal

ENV CRYPTOGRAPHY_DONT_BUILD_RUST=1

ARG password

ARG username

RUN echo $password# Pip password inputed: $password
RUN echo $username# Pip username inputed: $username

RUN yum install -y wget cmake nano python3-devel openssl-devel libffi-devel

RUN cd /usr/src && wget https://www.python.org/ftp/python/3.7.9/Python-3.7.9.tgz && tar xzf Python-3.7.9.tgz
RUN cd /usr/src && wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz && tar xzf Python-3.8.0.tgz

RUN cd /usr/src/Python-3.7.9 && ./configure --enable-optimizations && make altinstall
RUN cd /usr/src/Python-3.8.0 && ./configure --enable-optimizations && make altinstall

RUN python3   -m pip install --upgrade pip
RUN python3.7 -m pip install --upgrade pip
RUN python3.8 -m pip install --upgrade pip

RUN python3   -m pip install cryptography==3.3.2
RUN python3   -m pip install numpy scipy matplotlib cython pandas pyqt5 wheel vtk twine setuptools_rust
RUN python3.7 -m pip install numpy scipy matplotlib cython pandas pyqt5 wheel vtk twine setuptools_rust
RUN python3.8 -m pip install numpy scipy matplotlib cython pandas pyqt5 wheel vtk twine setuptools_rust


RUN echo Copying and cleaning files...
RUN mkdir GitProject && mkdir GitProject/PyMieSim
COPY . ./GitProject/PyMieSim/
RUN rm -rf /GitProject/PyMieSim/CMakeFiles /GitProject/PyMieSim/CMakeCache.txt
RUN mkdir /GitProject/PyMieSim/output
RUN rm -rf /GitProject/PyMieSim/dist/* /GitProject/PyMieSim/build/* /GitProject/PyMieSim/output/*
RUN rm /GitProject/PyMieSim/**/*.so -f
RUN rm /GitProject/PyMieSim/**/**/*.so -f


RUN echo Installing Boost dependencies...
RUN cp -r /GitProject/PyMieSim/extern/math/include/boost /usr/include


RUN echo Installing Fortran 90 dependencies...
RUN mkdir /GitProject/PyMieSim/extern/complex_bessel/build
RUN cd /GitProject/PyMieSim/extern/complex_bessel/build && cmake -DBUILD_TESTING=OFF .. && make install



RUN echo Compiling and uploading wheel for Python3.8
RUN cd /GitProject/PyMieSim && cmake . -DVERSION=3.8 && make clean && make all
RUN cd /GitProject/PyMieSim && python3.8 -m pip install -r requirements.txt
RUN cd /GitProject/PyMieSim/ && python3.8 setup.py install
RUN cd /GitProject/PyMieSim/ && python3.8 setup.py bdist_wheel
RUN auditwheel repair /GitProject/PyMieSim/dist/*38*.whl -w /GitProject/PyMieSim/output
RUN python3.8 -m twine upload --password $password --username $username --repository pypi /GitProject/PyMieSim/output/*38*manylinux2014*.whl
RUN rm -rf /GitProject/PyMieSim/output/*
RUN rm /GitProject/PyMieSim/CMakeCache.txt
RUN rm /GitProject/PyMieSim/**/*.so -f
RUN rm /GitProject/PyMieSim/**/**/*.so -f


RUN echo Compiling and uploading wheel for Python3.7
RUN cd /GitProject/PyMieSim && cmake . -DVERSION=3.7 && make clean && make all
RUN cd /GitProject/PyMieSim && python3.7 -m pip install -r requirements.txt
RUN cd /GitProject/PyMieSim/ && python3.7 setup.py install
RUN cd /GitProject/PyMieSim/ && python3.7 setup.py bdist_wheel
RUN auditwheel repair /GitProject/PyMieSim/dist/*37*.whl -w /GitProject/PyMieSim/output
RUN python3.8 -m twine upload --password $password --username $username --repository pypi /GitProject/PyMieSim/output/*37*manylinux2014*.whl
RUN rm -rf /GitProject/PyMieSim/output/*
RUN rm /GitProject/PyMieSim/CMakeCache.txt
RUN rm /GitProject/PyMieSim/**/*.so -f
RUN rm /GitProject/PyMieSim/**/**/*.so -f


RUN echo Compiling and uploading wheel for Python3.6
RUN cd /GitProject/PyMieSim && cmake . -DVERSION=3.6 && make clean && make all
RUN cd /GitProject/PyMieSim && python3.6 -m pip install -r requirements.txt
RUN cd /GitProject/PyMieSim/ && python3.6 setup.py install
RUN cd /GitProject/PyMieSim/ && python3.6 setup.py bdist_wheel
RUN auditwheel repair /GitProject/PyMieSim/dist/*36*.whl -w /GitProject/PyMieSim/output
RUN python3.8 -m twine upload --password $password --username $username --repository pypi /GitProject/PyMieSim/output/*36*manylinux2014*.whl
RUN rm -rf /GitProject/PyMieSim/output/*
RUN rm /GitProject/PyMieSim/CMakeCache.txt 
RUN rm /GitProject/PyMieSim/**/*.so -f
RUN rm /GitProject/PyMieSim/**/**/*.so -f
