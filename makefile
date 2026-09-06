.DEFAULT_GOAL := help

PYTHON ?= python3.13
BUILD_DIR ?= build
ROOT_DIR := $(CURDIR)
PYBIND11_DIR = $(shell $(PYTHON) -m pybind11 --cmakedir)
RELEASE_KIND := $(filter major minor patch,$(MAKECMDGOALS))

.PHONY: help configure build install quick rebuild editable clean quality test check release-check tag release major minor patch

help:
	@echo "PyMieSim development commands"
	@echo ""
	@echo "  make editable              Build and install an editable package"
	@echo "  make test                  Run the test suite"
	@echo "  make quality               Run static checks"
	@echo "  make check                 Run quality and tests"
	@echo "  make release-check         Check tag-derived release metadata"
	@echo "  make tag VERSION=vX.Y.Z    Create a release commit and annotated tag"
	@echo "  make release patch         Create and push the next patch release"
	@echo "  make release minor         Create and push the next minor release"
	@echo "  make release major         Create and push the next major release"

quality:
	$(PYTHON) -m ruff check PyMieSim tests
	$(PYTHON) -m mypy PyMieSim/gui/parsing.py PyMieSim/gui/schemas.py

test:
	$(PYTHON) -m pytest --config-file=pytest.ini

check: quality test

release-check:
	$(PYTHON) tools/check_release.py $(if $(VERSION),--version $(VERSION),)

tag:
	$(PYTHON) tools/release_tag.py "$(VERSION)"

release:
	@test "$(words $(RELEASE_KIND))" -eq 1 || { echo "usage: make release [patch|minor|major]" >&2; exit 2; }
	@set -eu; release_tag="$$($(PYTHON) tools/next_release_version.py $(RELEASE_KIND))"; \
	$(PYTHON) tools/release_tag.py "$$release_tag"; \
	git push origin HEAD "refs/tags/$$release_tag"

major minor patch:
	@:

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
