import os
import pytest
import sys

@pytest.fixture(params=['Fortran', 'Python'])
def backend_fixture(request):
    # Parametrizing the entire test suite: https://github.com/pytest-dev/pytest/issues/3196
    return request.param

@pytest.fixture(scope='function')
def pyprima_turn_on_debugging():
    # This is to force us to use DEBUGGING for a test, so that we get
    # code coverage. We definitely don't want to do this for the pycutest tests,
    # they are slow enough with USE_NAIVE_MATH set to True.
    from prima.backends.pyprima.common import consts
    consts.DEBUGGING[0] = True
    yield
    consts.DEBUGGING[0] = False
