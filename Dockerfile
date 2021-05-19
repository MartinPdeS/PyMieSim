FROM quay.io/pypa/manylinux2014_x86_64

ENV QT_QPA_PLATFORM minimal

ENV CRYPTOGRAPHY_DONT_BUILD_RUST=1

ARG password

ARG username

RUN yum install -y wget cmake nano python3-devel openssl-devel libffi-devel

RUN cd /usr/src && wget https://www.python.org/ftp/python/3.7.9/Python-3.7.9.tgz && tar xzf Python-3.7.9.tgz

RUN cd /usr/src && wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz && tar xzf Python-3.8.0.tgz

RUN cd /usr/src/Python-3.7.9 && ./configure --enable-optimizations && make altinstall

RUN cd /usr/src/Python-3.8.0 && ./configure --enable-optimizations && make altinstall

RUN ls -ls /usr/bin/python*

RUN mkdir GitProject && mkdir GitProject/PyMieSim

COPY . ./GitProject/PyMieSim/



RUN rm -rf dist/* build/*

RUN pip3 install -U pip

RUN python3.6 -m pip install --upgrade pip
RUN python3.7 -m pip install --upgrade pip
RUN python3.8 -m pip install --upgrade pip

RUN python3.6 -m pip install numpy scipy matplotlib cython pandas pyqt5 wheel vtk
RUN python3.7 -m pip install numpy scipy matplotlib cython pandas pyqt5 wheel vtk
RUN python3.8 -m pip install numpy scipy matplotlib cython pandas pyqt5 wheel vtk

RUN wget -c 'http://sourceforge.net/projects/boost/files/boost/1.75.0/boost_1_75_0.tar.bz2'

RUN tar xf boost_1_75_0.tar.bz2

RUN cp -r boost_1_75_0/boost /usr/include

RUN cd GitProject && git clone https://github.com/joeydumont/complex_bessel.git

RUN mkdir GitProject/complex_bessel/build

RUN cd GitProject/complex_bessel/build && cmake -DBUILD_TESTING=OFF .. && make install

RUN rm -rf GitProject/PyMieSim/CMakeFiles GitProject/PyMieSim/CMakeCache.txt

RUN cd GitProject/PyMieSim && cmake . && make all

RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"

RUN python3 get-pip.py

RUN cd GitProject/PyMieSim && pip3 install -r requirements.txt

RUN cd GitProject/PyMieSim && python3 setup.py install

RUN mkdir GitProject/PyMieSim/output

RUN python3 get-pip.py

RUN python3 -m pip install setuptools_rust wheel

RUN python3 -m pip install twine

RUN echo Pip password inputed: $password

RUN echo Pip username inputed: $username




RUN echo Compiling wheel for Python3.6

#RUN python3 -m pip wheel GitProject/PyMieSim -w GitProject/PyMieSim/output

RUN python3 GitProject/PyMieSim/setup.py bdist_wheel



RUN auditwheel repair GitProject/PyMieSim/dist/* -w GitProject/PyMieSim/output

RUN python3 -m twine upload --password $password --username $username --repository pypi GitProject/PyMieSim/output/PyMieSim*manylinux2014*




RUN echo Compiling wheel for Python3.7

RUN pip3.7 install wheel vtk

RUN python3.7 -m pip wheel GitProject/PyMieSim -w GitProject/PyMieSim/output

RUN auditwheel repair GitProject/PyMieSim/output/PyMieSim*cp37* -w GitProject/PyMieSim/output

RUN python3 -m twine upload --password $password --username $username --repository pypi GitProject/PyMieSim/output/PyMieSim*manylinux2014*



RUN echo Compiling wheel for Python3.8

RUN pip3.8 install wheel vtk

RUN python3.8 -m pip wheel /GitProject/PyMieSim -w /GitProject/PyMieSim/output

RUN auditwheel repair /GitProject/PyMieSim/output/PyMieSim*cp38* -w /GitProject/PyMieSim/output

RUN python3 -m twine upload --password $password --username $username --repository pypi /GitProject/PyMieSim/output/PyMieSim*manylinux2014*
