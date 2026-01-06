from prima import minimize, LinearConstraint as LC, NonlinearConstraint as NLC
import numpy as np
import pytest

def test_eliminate_fixed_bounds():
    # Test the logic for detecting and eliminating fixed bounds

    def f(x):
        return np.sum(x**2)

    lb = [-1, None, 1, None, -0.5]
    ub = [-0.5, -0.5, None, None, -0.5]
    bounds = [(a, b) for a, b in zip(lb, ub)]    
    res = minimize(f, x0=np.array([1, 2, 3, 4, 5]), bounds=bounds)
    assert np.allclose(res.x, np.array([-0.5, -0.5, 1, 0, -0.5]), atol=1e-3)
    assert np.allclose(res.fun, 1.75, atol=1e-3)


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
    assert np.isclose(res.fun, 51, atol=1e-6, rtol=1e-6)


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
    assert np.isclose(res.fun, 10, atol=1e-6, rtol=1e-6)


def test_infeasible_bounds():
    def f(x):
        return np.sum(x**2)

    lb = [1, None, None]
    ub = [-1, None, None]
    bounds = [(a, b) for a, b in zip(lb, ub)]
    with pytest.raises(AssertionError) as excinfo:
        minimize(f, x0=np.array([1, 2, 3]), bounds=bounds)
    assert str(excinfo.value) == "Some of the provided bounds are infeasible. infeasible=array([ True, False, False]) lb[infeasible]=array([1.]), ub[infeasible]=array([-1.])"


def test_infeasible_bounds_nan():
    def f(x):
        return np.sum(x**2)

    lb = [np.nan, None, None]
    ub = [-1, None, None]
    bounds = [(a, b) for a, b in zip(lb, ub)]
    with pytest.raises(AssertionError) as excinfo:
        minimize(f, x0=np.array([1, 2, 3]), bounds=bounds)
    assert str(excinfo.value) == "Some of the provided bounds are infeasible. infeasible=array([ True, False, False]) lb[infeasible]=array([nan]), ub[infeasible]=array([-1.])"


@pytest.mark.parametrize('with_Aineq', [False, True])
@pytest.mark.parametrize('with_Aeq', [False, True])
@pytest.mark.parametrize('with_nlcon', [False, True])
def test_all_fixed(with_nlcon, with_Aeq, with_Aineq):
    def f(x):
        return np.sum(x**2)

    lb = [-1, 2, 4]
    ub = [-1, 2, 4]
    bounds = [(a, b) for a, b in zip(lb, ub)]
    nlc = NLC(lambda x: x[2]**2, lb=-np.inf, ub=0)
    lc_eq = LC([0, 0, 1], lb=21, ub=21)
    lc_ineq = LC([0, 0, 1], lb=22, ub=23)
    constraints = []
    if with_nlcon:
        constraints.append(nlc)
    if with_Aeq:
        constraints.append(lc_eq)
    if with_Aineq:
        constraints.append(lc_ineq)

    result = minimize(f, x0=np.array([1, 2, 3]), bounds=bounds, constraints=constraints)
    assert all(result.x == [-1, 2, 4])
    assert result.fun == 21
    assert result.nfev == 1
    if not with_nlcon and not with_Aeq and not with_Aineq:
        assert np.array_equal(result.nlconstr, [])
        assert result.maxcv == 0
    if not with_nlcon and not with_Aeq and with_Aineq:
        assert np.array_equal(result.nlconstr, [])
        assert result.maxcv == 18
    if not with_nlcon and with_Aeq and not with_Aineq:
        assert np.array_equal(result.nlconstr, [])
        assert result.maxcv == 17
    if not with_nlcon and with_Aeq and with_Aineq:
        assert np.array_equal(result.nlconstr, [])
        assert result.maxcv == 18
    if with_nlcon and not with_Aeq and not with_Aineq:
        assert np.array_equal(result.nlconstr, [16])
        assert result.maxcv == 16
    if with_nlcon and not with_Aeq and with_Aineq:
        assert np.array_equal(result.nlconstr, [16])
        assert result.maxcv == 18
    if with_nlcon and with_Aeq and not with_Aineq:
        assert np.array_equal(result.nlconstr, [16])
        assert result.maxcv == 17
    if with_nlcon and with_Aeq and with_Aineq:
        assert np.array_equal(result.nlconstr, [16])
        assert result.maxcv == 18
