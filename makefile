PYTHON ?= python3.13
BUILD_DIR ?= build
ROOT_DIR := $(CURDIR)
PYBIND11_DIR := $(shell $(PYTHON) -m pybind11 --cmakedir)

.PHONY: configure build install quick rebuild editable clean

configure:
	cmake -S . -B $(BUILD_DIR) \
		-Dpybind11_DIR="$(PYBIND11_DIR)" \
		-DPython_EXECUTABLE="$$(which $(PYTHON))" \
		-DCMAKE_INSTALL_PREFIX="$(ROOT_DIR)"

build:
	cmake --build $(BUILD_DIR) -j

install:
	cmake --install $(BUILD_DIR)

uninstall:
	$(PYTHON) -m pip uninstall -y PyMieSim

quick: configure build install

rebuild: configure build install

editable:
	$(PYTHON) -m pip install --no-build-isolation -Cbuild-dir=build -Ceditable.rebuild=false -Ceditable.mode=inplace -e .

clean:
	rm -rf $(BUILD_DIR)