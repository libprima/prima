# On some platforms in CI we are not able to install scipy, and in that
# case we should skip this test. Note that pdfo depends on scipy.
import pytest
from packaging import version
from test_combining_constraints import test_providing_bounds_and_linear_and_nonlinear_constraints


def test_prima():
    from prima import minimize, NonlinearConstraint as NLC, LinearConstraint as LC, Bounds
    test_providing_bounds_and_linear_and_nonlinear_constraints(minimize, NLC, LC, Bounds)


# Despite the fact that we are using the pdfo function, we still get this warning because the pdfo
# function itself calls the cobyla function. For now we suppress this warning.
@pytest.mark.filterwarnings("ignore:The `cobyla` function is deprecated. Use the `pdfo` function")
def test_pdfo():
    pdfo = pytest.importorskip("pdfo")
    if version.parse(pdfo.__version__) < version.parse("2.0.0"):
        pytest.skip("pdfo version too old for this test (it does not accept radius_init and radius_final as options)")
    scipy = pytest.importorskip("scipy")
    if version.parse(scipy.__version__) < version.parse("1.11.0"):
        pytest.skip("scipy version too old for this test (its version of COBYLA does not accept bounds)")
    numpy = pytest.importorskip("numpy")
    if version.parse(numpy.__version__) >= version.parse("2.0.0"):
        pytest.skip("numpy version too new for this test (pdfo does not yet support numpy v2)")

    from pdfo import pdfo
    from scipy.optimize import NonlinearConstraint as NLC, LinearConstraint as LC, Bounds
    test_providing_bounds_and_linear_and_nonlinear_constraints(pdfo, NLC, LC, Bounds, package='pdfo')
