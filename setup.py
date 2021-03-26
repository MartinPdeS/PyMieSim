#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev





# command: sudo python setup.py upload_sphinx
import io
import os
import sys
from shutil import rmtree
from setuptools import setup, Extension
#from Cython.Distutils import build_ext

from setuptools import find_packages, setup, Command
from setuptools.command.build_ext import build_ext as _build_ext


# Package meta-data.
NAME            = 'PyMieSim'
DESCRIPTION     = 'A package for light scattering simulations.'
URL             = 'https://github.com/MartinPdeS/PyMieSim'
EMAIL           = 'Martin.poinsinet.de.sivry@gmail.com'
AUTHOR          = 'Martin Poinsinet de Sivry',
REQUIRES_PYTHON = '>=3.6.0'
VERSION         = '0.1.11'

# What packages are required for this module to be executed?
REQUIRED = ['scipy',
            'matplotlib',
            'pandas',
            'mayavi']


class get_pybind11_include(object):
    """Defer numpy.get_include() until after numpy is installed."""

    def __str__(self):
        import pybind11
        return pybind11.get_include()

class get_numpy_include(object):
    """Defer numpy.get_include() until after numpy is installed."""

    def __str__(self):
        import numpy
        return numpy.get_include()



# What packages are optional?
EXTRAS = {
    # 'fancy feature': ['django'],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))


macro = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]
compile_args=["-std=c++14",
              '-fopenmp',
              '-lboost_filesystem',
              '-lboost_system',
              '-O3',
              '-march=native']

link_args=["-std=c++14",
           '-fopenmp',
           '-lboost_filesystem',
           '-lboost_system',
           '-O3',
           '-march=native']

ext_modules = [
                Extension(name               = "PyMieSim._Coupling",
                          sources            = ["PyMieSim/Coupling.pyx"],
                          include_dirs        = [get_numpy_include()],
                          language            = "c++",
                          define_macros       = macro,
                          extra_compile_args  = compile_args,
                          extra_link_args     = link_args),

                Extension(name               = "PyMieSim.LMT.Scatterer",
                          sources            = ["PyMieSim/LMT/cpp/interface.cpp"],
                          include_dirs       = [get_numpy_include(), get_pybind11_include()],
                          extra_compile_args = compile_args,
                          language           = "c++"),

                Extension(name                = 'PyMieSim.GLMT.Scatterer',
                          sources             = ['PyMieSim/GLMT/cpp/Interface.cpp'],
                          include_dirs        = [get_numpy_include(), get_pybind11_include()],
                          extra_compile_args  = compile_args,
                          language            = 'c++'),

                Extension(name         = 'PyMieSim.Fibonacci',
                          sources      = ['PyMieSim/FibonnaciMesh.cpp'],
                          include_dirs = [get_numpy_include(), get_pybind11_include()],
                          extra_compile_args  = compile_args,
                          language     = 'c++'),

                Extension(name         = 'PyMieSim.GLMT.GaussianBeam',
                          sources      = ['PyMieSim/GLMT/cpp/GaussianBeam.cpp'],
                          include_dirs = [get_numpy_include(), get_pybind11_include()],
                          language     = 'c++',
                          extra_compile_args  = compile_args)
                ]


# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()


# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    setup_requires=['numpy', 'pybind11','cython'],
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],

    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    dependency_links = ["fibermodes @ git+https://github.com/cbrunet/fibermodes.git#egg=fibermodes-0.2.0"],
    include_package_data = True,
    ext_modules          = ext_modules,
    license='MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ],
    # $ setup.py publish support.
    cmdclass={
        'upload': UploadCommand,
        #'build_ext': build_ext
    },
)
