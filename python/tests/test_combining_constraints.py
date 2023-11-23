from prima import minimize as prima_minimize, NonlinearConstraint as prima_NLC, LinearConstraint as prima_LC, Bounds as prima_Bounds
import numpy as np
from objective import fun


def test_providing_linear_and_nonlinear_constraints(capfd):
    nlc = prima_NLC(lambda x: x[0]**2, lb=[25], ub=[100])
    lc = prima_LC(np.array([1,1]), lb=10, ub=15)
    x0 = [0, 0]
    res = prima_minimize(fun, x0, constraints=[nlc, lc])
    assert np.isclose(res.x[0], 5.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 4.5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 0.5, atol=1e-6, rtol=1e-6)
    outerr = capfd.readouterr()
    assert outerr.out == "Nonlinear constraints detected, applying COBYLA\n"
    assert outerr.err == ''


def test_providing_bounds_and_linear_constraints(capfd):
    lc = prima_LC(np.array([1,1]), lb=10, ub=15)
    bounds = prima_Bounds(1, 1)
    x0 = [0, 0]
    res = prima_minimize(fun, x0, constraints=lc, bounds=bounds)
    assert np.isclose(res.x[0], 1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 9, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 41, atol=1e-6, rtol=1e-6)
    outerr = capfd.readouterr()
    assert outerr.out == "Linear constraints detected without nonlinear constraints, applying LINCOA\n"
    assert outerr.err == ''


def test_providing_bounds_and_nonlinear_constraints(capfd):
    nlc = prima_NLC(lambda x: x[0]**2, lb=[25], ub=[100])
    bounds = prima_Bounds([None, 1], [None, 1])
    x0 = [0, 0]
    res = prima_minimize(fun, x0, constraints=nlc, bounds=bounds)
    assert np.isclose(res.x[0], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.fun, 9, atol=1e-6, rtol=1e-6)
    outerr = capfd.readouterr()
    assert outerr.out == "Nonlinear constraints detected, applying COBYLA\n"
    assert outerr.err == ''


# This test is re-used for the compatibility tests, hence the extra arguments and their
# default values
def test_providing_bounds_and_linear_and_nonlinear_constraints(capfd, minimize=prima_minimize, NLC=prima_NLC, LC=prima_LC, Bounds=prima_Bounds, package='prima'):
    # This test needs a 3 variable objective function so that we can check that the
    # bounds and constraints are all active
    def newfun(x):
        return fun(x[0:2]) + (x[2] - 3)**2
    nlc = NLC(lambda x: x[0]**2, lb=[25], ub=[100])
    bounds = Bounds([-np.inf, 1, -np.inf], [np.inf, 1, np.inf])
    lc = LC(np.array([1,1,1]), lb=10, ub=15)
    x0 = [0, 0, 0]
    # macOS seems to stop just short of the optimal solution, so we help it along by
    # taking a larger initial trust region radius and requiring a smaller final radius
    # before stopping. The different packages have different names for these options.
    if package == 'prima':
        options = {'rhobeg': 2, 'rhoend': 1e-8}
        method = None
    elif package == 'pdfo':
        options = {'radius_init': 2, 'radius_final': 1e-8}
        method = None
    elif package == 'scipy':
        options = {'rhobeg': 2, 'tol': 1e-8}
        # PDFO and PRIMA will select COBYLA but SciPy may select something else, so we tell it to select COBYLA
        method = 'COBYLA'
    else:
        # Since this is test infrastructure under the control of the developers we
        # should never get here except for a typo or something like that
        raise ValueError(f"Unknown package: {package}")

    res = minimize(newfun, x0, method=method, constraints=[nlc, lc], bounds=bounds, options=options)
    
    # 32 bit builds of PRIMA reach the optimal solution with the same level of precision as 64 bit builds
    # so we lower the atol/rtol to 1e-3 so that 32 bit builds will pass.
    assert np.isclose(res.x[0], 5.5, atol=1e-3, rtol=1e-3)
    assert np.isclose(res.x[1], 1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[2], 3.5, atol=1e-3, rtol=1e-3)
    assert np.isclose(res.fun, 9.5, atol=1e-3, rtol=1e-3)
    if package == 'prima':
        outerr = capfd.readouterr()
        assert outerr.out == "Nonlinear constraints detected, applying COBYLA\n"
        assert outerr.err == ''
