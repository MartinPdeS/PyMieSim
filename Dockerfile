FROM quay.io/pypa/manylinux2014_x86_64

ENV QT_QPA_PLATFORM minimal

ENV CRYPTOGRAPHY_DONT_BUILD_RUST=1

ARG password

ARG username

RUN mkdir GitProject && mkdir GitProject/PyMieSim

COPY . ./GitProject/PyMieSim/



RUN yum install -y wget cmake nano python3-devel openssl-devel

RUN pip3 install -U pip

RUN pip3 install numpy scipy matplotlib cython pandas pyqt5

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


RUN python3 -m pip wheel GitProject/PyMieSim -w GitProject/PyMieSim/output

RUN auditwheel repair GitProject/PyMieSim/output/PyMieSim*whl -w GitProject/PyMieSim/output

RUN python3 -m twine upload --password $password --username $username --repository pypi GitProject/PyMieSim/output/PyMieSim*manylinux2014*
