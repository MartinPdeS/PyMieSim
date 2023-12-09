# simple makefile to simplify repetitive build env management tasks under posix

# caution: testing won't work on windows, see README

PYTHON ?= python
PYTEST ?= pytest

test:
	rm -rf coverage .coverage
	$(PYTEST) --showlocals -v --cov=PyMieSim --cov-report=html:coverage

trailing-spaces:
	find PyMieSim -name "*.py" -exec perl -pi -e 's/[ \t]*$$//' {} \;

doc:
	$(MAKE) -C docs html

doc-noplot:
	$(MAKE) -C docs html-noplot

clean:
	rm -rf build PyMieSim.egg-info dist .pytest_cache PyMieSim/binary/*.a PyMieSim/binary/*.so *.whl wheel_house .eggs