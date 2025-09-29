from pyprima import minimize, LinearConstraint as LC, NonlinearConstraint as NLC
import numpy as np

def test_eliminate_fixed_bounds():
    # Test the logic for detecting and eliminating fixed bounds

    def f(x):
        return np.sum(x**2)

    lb = [-1, None, 1, None, -0.5]
    ub = [-0.5, -0.5, None, None, -0.5]
    bounds = [(a, b) for a, b in zip(lb, ub)]    
    res = minimize(f, x0=np.array([1, 2, 3, 4, 5]), bounds=bounds)
    assert np.allclose(res.x, np.array([-0.5, -0.5, 1, 0, -0.5]), atol=1e-3)
    assert np.allclose(res.f, 1.75, atol=1e-3)


def test_eliminate_fixed_bounds_with_linear_constraints():
    # Ensure that the logic for fixed bounds also modifies linear constraints
    # appropriately

    def f(x):
        return np.sum(x**2)

    lb = [-1, None, None]
    ub = [-1, None, None]
    bounds = [(a, b) for a, b in zip(lb, ub)]
    # Internally, the algorithm should modify lc to be a 1x2 matrix instead of 1x3,
    # since it will modify x0 and the objective function to eliminate the first
    # variable. If it fails to modify lc, we will get shape mismatch errors.
    lc = LC(np.array([1, 1, 1]), lb=9, ub=15)
    res = minimize(f, x0=np.array([1, 1, 1]), constraints=lc, bounds=bounds)
    assert np.isclose(res.x[0], -1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[2], 5, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.f, 51, atol=1e-6, rtol=1e-6)


def test_eliminate_fixed_bounds_with_nonlinear_constraints():
    # Ensure that the logic for fixed bounds also modifies the nonlinear constraint
    # function appropriately

    def f(x):
        return np.sum(x**2)

    lb = [-1, None, None]
    ub = [-1, None, None]
    bounds = [(a, b) for a, b in zip(lb, ub)]
    x0 = np.array([1, 1, 1])
    # Have the nonlinear constraint function operate on the last element of x, but be
    # explicit about the length of x. This ensures that the test is still valid if the
    # fixed bound is removed. If we simply used x[-1] this test would pass but it
    # wouldn't actually test if we had properly modified the nonlinear constraint
    # function after removing the fixed bounds
    nlc = NLC(lambda x: x[len(x0)-1]**2, lb=9, ub=15)
    res = minimize(f, x0=x0, constraints=nlc, bounds=bounds)
    assert np.isclose(res.x[0], -1, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[1], 0, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.x[2], 3, atol=1e-6, rtol=1e-6)
    assert np.isclose(res.f, 10, atol=1e-6, rtol=1e-6)