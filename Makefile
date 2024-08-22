# simple makefile to simplify repetitive build env management tasks under posix

# caution: testing won't work on windows, see README

PYTHON ?= python
PYTEST ?= pytest

test:
	rm -rf coverage .coverage
	$(PYTEST) --showlocals -v --cov=PyMieSim --cov-report=html:coverage tests

trailing-spaces:
	find PyMieSim -name "*.py" -exec perl -pi -e 's/[ \t]*$$//' {} \;

doc:
	$(MAKE) -C docs html

doc-noplot:
	$(MAKE) -C docs html-noplot

clean_pycache:
	find . -name __pycache__ -prune -exec rm -rf {} \;

clean_shared_object:
	find . -name "*.so" -prune -exec rm -rf {} \;

clean: clean_pycache clean_shared_object
	rm -rf build *.egg-info dist .pytest_cache *.whl wheel_house .eggs .coverage