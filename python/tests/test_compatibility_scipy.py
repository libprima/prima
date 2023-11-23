# On some platforms in CI we are not able to install scipy, and in that
# case we should skip this test. Note that pdfo depends on scipy.
import pytest
from packaging import version

from test_combining_constraints import test_providing_bounds_and_linear_and_nonlinear_constraints

def test_prima(capfd):
    from prima import minimize, NonlinearConstraint as NLC, LinearConstraint as LC, Bounds
    test_providing_bounds_and_linear_and_nonlinear_constraints(capfd, minimize, NLC, LC, Bounds)


def test_scipy():
    scipy = pytest.importorskip("scipy")
    if version.parse(scipy.__version__) < version.parse("1.11.0"):
        pytest.skip("scipy version too old for this test (its version of COBYLA does not accept bounds)")
    from scipy.optimize import minimize, NonlinearConstraint as NLC, LinearConstraint as LC, Bounds
    test_providing_bounds_and_linear_and_nonlinear_constraints(None, minimize, NLC, LC, Bounds, package="scipy")
