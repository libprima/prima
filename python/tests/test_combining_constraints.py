from prima import minimize as prima_minimize, NonlinearConstraint as prima_NLC, LinearConstraint as prima_LC, Bounds as prima_Bounds
import numpy as np
from objective import fun


def test_providing_linear_and_nonlinear_constraints():
    nlc = prima_NLC(lambda x: x[0]**2, lb=[25], ub=[100])
    lc = prima_LC(np.array([1,1]), lb=10, ub=15)
    x0 = [0, 0]
    res = prima_minimize(fun, x0, constraints=[nlc, lc])
    assert np.isclose(res.x[0], 5.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 4.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 0.5, atol=1e-6, rtol=1e-6)
    assert res.method == "cobyla"


def test_providing_bounds_and_linear_constraints():
    lc = prima_LC(np.array([1,1]), lb=10, ub=15)
    bounds = prima_Bounds(1, 1)
    x0 = [0, 0]
    res = prima_minimize(fun, x0, constraints=lc, bounds=bounds)
    assert np.isclose(res.x[0], 1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 9, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 41, atol=1e-6, rtol=1e-6)
    assert res.method == "lincoa"


def test_providing_bounds_and_nonlinear_constraints():
    nlc = prima_NLC(lambda x: x[0]**2, lb=[25], ub=[100])
    bounds = prima_Bounds([None, 1], [None, 1])
    x0 = [6, 1]  # Unfortunately the test is very fragile if we do not start near the optimal point
    res = prima_minimize(fun, x0, constraints=nlc, bounds=bounds)
    assert np.isclose(res.x[0], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 9, atol=1e-6, rtol=1e-6)
    assert res.method == "cobyla"


# This test is re-used for the compatibility tests, hence the extra arguments and their
# default values
def test_providing_bounds_and_linear_and_nonlinear_constraints(minimize=prima_minimize, NLC=prima_NLC, LC=prima_LC, Bounds=prima_Bounds, package='prima'):
    # This test needs a 3 variable objective function so that we can check that the
    # bounds and constraints are all active
    def newfun(x):
        return fun(x[0:2]) + (x[2] - 3)**2
    nlc = NLC(lambda x: x[0]**2, lb=[25], ub=[100])
    bounds = Bounds([-np.inf, 1, -np.inf], [np.inf, 1, np.inf])
    lc = LC(np.array([1,1,1]), lb=10, ub=15)
    x0 = [6, 1, 3.5]  # The test becomes very fragile if we do not start near the optimal point
    # PDFO and PRIMA will select COBYLA but SciPy may select something else, so we tell it to select COBYLA
    method = 'COBYLA' if package == 'scipy' else None
    res = minimize(newfun, x0, method=method, constraints=[nlc, lc], bounds=bounds)

    # 32 bit builds of PRIMA reach the optimal solution with the same level of precision as 64 bit builds
    # so we lower the atol/rtol to 1e-3 so that 32 bit builds will pass.
    assert np.isclose(res.x[0], 5.5, atol=1e-3, rtol=1e-3)
    assert np.isclose(res.x[1], 1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[2], 3.5, atol=1e-3, rtol=1e-3)
    assert np.isclose(res.fun, 9.5, atol=1e-3, rtol=1e-3)
    if package == 'prima' or package == 'pdfo':
        assert res.method == "cobyla"
