.PHONY: all test test-no-cmake coverage-html

all: coverage-html

test:
	SKBUILD_CMAKE_BUILD_TYPE=Debug uv run --reinstall pytest -s --cov=. --cov-branch --durations=20 tests

test-no-cmake:
	SKBUILD_WHEEL_CMAKE=0 COVERAGE_FILE=.coverage.2 uv run --reinstall pytest --cov=. --cov-branch tests/test_without_pybindings.py

coverage-html: test test-no-cmake
	uvx coverage combine --append
	uvx coverage html -d $$(git rev-parse --show-toplevel)/prima_htmlcov

test-quick:
	uv run pytest -s --cov=. --cov-branch --durations=20 --ignore=tests/pyprima/test_pycutest.py tests

coverage-quick: test-quick
	uvx coverage html -d $$(git rev-parse --show-toplevel)/prima_htmlcov