import pytest
# This exists mainly for those CI tests in which cutest/pycutest are not installed.
import os
import sys
pytestmark = [
    pytest.mark.skipif(os.getenv('SKBUILD_CMAKE_BUILD_TYPE') != "Debug",
        reason="CUTEST tests must be run against Fortran compiled in Debug mode"),
    pytest.mark.skipif(sys.platform not in ['linux', 'darwin'],
        reason="pycutest is not supported on Windows")
]

from .load_cutest_problem import load_cutest_problem
import prima
from prima import minimize, Bounds, LinearConstraint, NonlinearConstraint
import numpy as np


'''
This module tests various problem from the CUTEST set in order to really stress test
the implementation and also cover some code not covered by the naive tests in
test_end_to_end.py. The list is semi-arbitrary, some of these helped to find bugs
when testing the Python implementation against the Fortran one.

HS103 and MGH10LS provide coverage for the code in trustregion.py which scales the
problem if A contains large values.

HS102 is a very important problem. After the performance profiles were nearly identical
this problem led to non-identical results between Python and Fortran. The culprit was
found in getcpen - input values were being modified because they were being passed by
reference. Surprisingly this did not lead to issues in 99% of cases. After this
discovery the input values were copied being allowing getcpen to continue.

POLAK3, with the given options, provides coverage for the code to take an inverse of a
nontriangular matrix (via updatexfc)
'''


@pytest.fixture(autouse=True, scope='module')
def set_comparing():
    # This is a hack to force these tests to use manual math instead of optimized
    # numpy or other routines. This will ensure we get the same results as compared
    # to Fortran when compiled in debug mode.
    prima.pyprima.common.linalg.USE_NAIVE_MATH = True
    yield
    prima.pyprima.common.linalg.USE_NAIVE_MATH = False


def get_constraints(constraints_in):
    constraints_out = []
    if constraints_in['b_ub'].size > 0:
        constraints_out.append(LinearConstraint(constraints_in['a_ub'], -np.inf, constraints_in['b_ub']))
    if constraints_in['b_eq'].size > 0:
        constraints_out.append(LinearConstraint(constraints_in['a_eq'], constraints_in['b_eq'], constraints_in['b_eq']))
    if constraints_in['m_nonlinear_ub'] > 0:
        constraints_out.append(NonlinearConstraint(constraints_in['c_ub'], -np.inf, np.zeros(constraints_in['m_nonlinear_ub'])))
    if constraints_in['m_nonlinear_eq'] > 0:
        constraints_out.append(NonlinearConstraint(constraints_in['c_eq'], np.zeros(constraints_in['m_nonlinear_eq']), np.zeros(constraints_in['m_nonlinear_eq'])))
    return constraints_out


@pytest.mark.parametrize("name", [#'TENBARS1', 'ERRINBAR',
                                  'PALMER2C', 'PALMER3B',
                                  'HS103', 'CRESC4', 'MGH10LS', 'TFI3',
                                  'BIGGS3', 'BIGGS6', 'DEGENLPB', 'HS102',
                                  'MISRA1ALS', 'POLAK3'])
def test_cutest(name):
    fun, x0, lb, ub, constraints_in = load_cutest_problem(name)
    constraints = get_constraints(constraints_in) if constraints_in is not None else None
    bounds = Bounds(lb, ub)
    options = {'backend': 'Python'}
    if name == "POLAK3":
        options['maxfev'] = 271*len(x0)
        options['rhobeg'] = 2.71828
        options['rhoend'] = 3.14e-4
    options_copy = options.copy()
    python_result = minimize(fun, x0, method='cobyla', constraints=constraints, bounds=bounds, options=options)
    options_copy['backend'] = 'Fortran'
    fortran_result = minimize(fun, x0, method='cobyla', constraints=constraints, bounds=bounds, options=options_copy)
    assert np.allclose(python_result.x, fortran_result.x, atol=1e-15)
    assert np.isclose(python_result.fun, fortran_result.fun, atol=1e-15)
    assert np.allclose(python_result.nlconstr, fortran_result.nlconstr, atol=1e-15)
    assert python_result.nfev == fortran_result.nfev
