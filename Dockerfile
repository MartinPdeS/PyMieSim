FROM quay.io/pypa/manylinux2014_x86_64

ARG password

ARG username

RUN export 

RUN echo Pip password inputed: $password

RUN echo Pip username inputed: $username

RUN yum install -y wget cmake nano python3-devel

RUN wget -c 'http://sourceforge.net/projects/boost/files/boost/1.75.0/boost_1_75_0.tar.bz2'

RUN tar xf boost_1_75_0.tar.bz2

RUN cp -r boost_1_75_0/boost /usr/include

RUN mkdir GitProject

RUN cd GitProject && git clone https://github.com/joeydumont/complex_bessel.git 

RUN mkdir GitProject/complex_bessel/build

RUN cd GitProject/complex_bessel/build && cmake -DBUILD_TESTING=OFF .. && make install

RUN cd GitProject && git clone https://github.com/MartinPdeS/PyMieSim.git

RUN rm -rf GitProject/PyMieSim/CMakeFiles GitProject/PyMieSim/CMakeCache.txt 

RUN cd GitProject/PyMieSim && cmake . && make all 

RUN mkdir GitProject/PyMieSim/output

RUN /opt/python/cp37-cp37m/bin/pip install numpy vtk twine

RUN /opt/python/cp37-cp37m/bin/pip install mayavi

RUN /opt/python/cp37-cp37m/bin/pip wheel GitProject/PyMieSim -w GitProject/PyMieSim/output

RUN auditwheel repair GitProject/PyMieSim/output/PyMieSim*whl -w /output

RUN /opt/python/cp37-cp37m/bin/python3 -m twine upload --password $password --username $username --repository pypi GitProject/PyMieSim/dist/*



